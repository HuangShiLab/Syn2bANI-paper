# Reproducible environment for the Syn2bANI paper analysis and figure pipeline.
# The Syn2bANI binary itself is Rust and is built from
# https://github.com/HuangShiLab/Syn2bANI (not required for re-plotting).
#
# Build:  docker build -t syn2bani-paper .
# Run:    docker run --rm -v $PWD:/work -w /work syn2bani-paper \
#           python3 scripts/generate_fig8_cagpai.py

FROM continuumio/miniconda3:24.7.1-0

WORKDIR /opt
COPY environment.yml /opt/environment.yml
RUN conda env create -f /opt/environment.yml && conda clean -afy

ENV PATH=/opt/conda/envs/syn2bani-paper/bin:$PATH
WORKDIR /work

CMD ["bash"]
