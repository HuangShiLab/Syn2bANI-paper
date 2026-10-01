#!/usr/bin/env bash
# Lightweight status monitor for the local OAPGC census and HPC GTDB ANI pass.
# It records snapshots every 10 minutes and posts a macOS notification on
# completion, stall, or monitor shutdown conditions.
set -u

ROOT=/Users/macstudio/Downloads/OAPGC
WORK="$ROOT/census"
LOG="$ROOT/monitor_status.log"
EVENT="$ROOT/monitor_events.log"
INTERVAL="${INTERVAL:-600}"
LOCAL_EXPECTED=2120
GTDB_STRUCTURAL_EXPECTED=25582
GTDB_ANI_EXPECTED=22535

mkdir -p "$ROOT"
touch "$LOG" "$EVENT"

notify() {
  local msg="$1"
  printv "$(date -u +%FT%TZ) EVENT $msg" >> "$EVENT"
}

printv() { echo "$@" | tee -a "$LOG"; }

local_done=0
local_out=0
local_active=0
local_screen=0
local_free=0
struct_out=0
struct_ckpt=0
ani_out=0
ani_tmp=0
hpc_jobs='none'
hpc_state='UNKNOWN'
unchanged_ani=0
last_ani=-1
last_local=-1
local_complete_notified=0
ani_complete_notified=0
struct_complete_notified=0

while true; do
  ts=$(date -u +%FT%TZ)
  [ -f "$WORK/tasks.jsonl" ] && local_expected_now=$(wc -l < "$WORK/tasks.jsonl" | tr -d ' ') || local_expected_now="$LOCAL_EXPECTED"
  local_done=$(find "$WORK/checkpoints" -name '*.done' 2>/dev/null | wc -l | tr -d ' ')
  local_out=$(find "$WORK/outputs" -type f -name '*.tsv.gz' 2>/dev/null | wc -l | tr -d ' ')
  local_active=$(pgrep -f '/Users/macstudio/Downloads/Syn2b/target/release/syn2b synteny' 2>/dev/null | wc -l | tr -d ' ')
  screen -ls 2>/dev/null | grep -q '[.]oapgc' && local_screen=1 || local_screen=0
  local_free=$(df -g /Users/macstudio/Downloads 2>/dev/null | awk 'NR==2{print $4}')

  if (( local_done >= local_expected_now )); then
    local_state=COMPLETE
  elif (( local_active > 0 )); then
    local_state=RUNNING
  elif (( local_screen )); then
    local_state=RUNNING_WAITING
  else
    local_state=STALLED
  fi

  hpc_blob=$(ssh -o ConnectTimeout=20 hpc "
    W=/lustre1/g/aos_shihuang/data/syn2b_census;
    printf 'struct_out %s\n' \"\$(find \"\$W/outputs\" -maxdepth 1 -type f -name '*.tsv.gz' 2>/dev/null | wc -l)\";
    printf 'struct_ckpt %s\n' \"\$(find \"\$W/checkpoints\" -maxdepth 1 -type f -name '*.done' 2>/dev/null | wc -l)\";
    printf 'ani_out %s\n' \"\$(find \"\$W/ani\" -maxdepth 1 -type f -name '*.ani.tsv.gz' 2>/dev/null | wc -l)\";
    printf 'ani_tmp %s\n' \"\$(find \"\$W/ani\" -maxdepth 1 -type f -name '*.tmp' 2>/dev/null | wc -l)\";
    squeue -u \$USER -h -o '%i|%j|%T|%M|%L|%N' 2>/dev/null | grep census | head -10
  " 2>/dev/null)

  [ -z "$hpc_blob" ] && hpc_blob='struct_out 0
struct_ckpt 0
ani_out 0
ani_tmp 0
SSH_ERROR'
  struct_out=$(awk '$1=="struct_out"{print $2}' <<<"$hpc_blob")
  struct_ckpt=$(awk '$1=="struct_ckpt"{print $2}' <<<"$hpc_blob")
  ani_out=$(awk '$1=="ani_out"{print $2}' <<<"$hpc_blob")
  ani_tmp=$(awk '$1=="ani_tmp"{print $2}' <<<"$hpc_blob")
  hpc_jobs=$(grep '|' <<<"$hpc_blob" | head -5 || true)
  if grep -q 'census-ani|RUNNING' <<<"$hpc_jobs"; then hpc_state=RUNNING
  elif grep -q 'census-ani|PENDING' <<<"$hpc_jobs"; then hpc_state=PENDING
  elif grep -q SSH_ERROR <<<"$hpc_blob"; then hpc_state=SSH_ERROR
  else hpc_state=NOT_RUNNING; fi

  if (( struct_out >= GTDB_STRUCTURAL_EXPECTED && struct_ckpt >= GTDB_STRUCTURAL_EXPECTED )); then
    struct_state=COMPLETE
  else
    struct_state=INCOMPLETE
  fi
  if (( ani_out >= GTDB_ANI_EXPECTED )); then ani_state=COMPLETE; else ani_state=INCOMPLETE; fi

  printv "$ts LOCAL state=$local_state done=$local_done/$local_expected_now outputs=$local_out active_syn2b=$local_active screen=$local_screen free_GB=$local_free"
  printv "$ts HPC structural=$struct_state outputs=$struct_out checkpoints=$struct_ckpt ANI=$ani_state outputs=$ani_out/$GTDB_ANI_EXPECTED tmp=$ani_tmp state=$hpc_state jobs=${hpc_jobs//$'\n'/ ; }"

  if (( struct_out >= GTDB_STRUCTURAL_EXPECTED && struct_ckpt >= GTDB_STRUCTURAL_EXPECTED && ! struct_complete_notified )); then
    notify "GTDB structural census complete: $struct_out tasks"
    struct_complete_notified=1
  fi
  if (( local_done >= local_expected_now && ! local_complete_notified )); then
    notify "OAPGC structural census complete: $local_done/$local_expected_now tasks"
    local_complete_notified=1
  fi
  if (( ani_out >= GTDB_ANI_EXPECTED && ! ani_complete_notified )); then
    notify "GTDB ANI pass complete: $ani_out/$GTDB_ANI_EXPECTED clusters"
    ani_complete_notified=1
  fi

  if (( last_local >= 0 && last_local == local_done )) && [[ "$local_state" != COMPLETE ]]; then
    notify "Warning: OAPGC progress unchanged at $local_done/$local_expected_now"
  fi
  if (( last_ani >= 0 && last_ani == ani_out )) && [[ "$ani_state" != COMPLETE && "$hpc_state" == RUNNING ]]; then
    ((unchanged_ani++))
    if (( unchanged_ani == 6 )); then
      notify "Warning: GTDB ANI unchanged for ~1h at $ani_out/$GTDB_ANI_EXPECTED"
    fi
  else
    unchanged_ani=0
  fi
  last_local=$local_done
  last_ani=$ani_out

  if [[ "$local_state" == COMPLETE && "$struct_state" == COMPLETE && "$ani_state" == COMPLETE ]]; then
    notify "All monitored census tasks complete"
    break
  fi
  sleep "$INTERVAL"
done
