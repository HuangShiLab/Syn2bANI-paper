# FracMinHash scale sweep validation

## Accuracy by scale

| scale | n | mae | rmse | pearson_r | mean_shared_tags | median_shared_tags |
| --- | --- | --- | --- | --- | --- | --- |
| four | 43312 | 0.10946252640245933 | 0.19587819449184546 | 0.17705929456187144 | 561.2597201699298 | 314.0 |
| fmh750 | 43329 | 0.11146664290762955 | 0.19678408915923598 | 0.17595570058909088 | 422.9030903090309 | 254.0 |
| fmh1582 | 43251 | 0.11801376958736301 | 0.20061302613680065 | 0.16956867266644085 | 200.58192874153198 | 120.0 |

## Error-variance model

Fitted `Var(residual) = c0 + c1 * p*(1-p)/m`, where `p` is dnadiff inverted fraction and `m` is shared_tags.

| scale | c0 | c1 | r2 | n |
| --- | --- | --- | --- | --- |
| four | 0.03979253845817865 | -0.8362051900836897 | 0.0004879331512321894 | 43312 |
| fmh750 | 0.04328458144593508 | -2.959839972753355 | 0.0026200842494830123 | 43329 |
| fmh1582 | 0.04224754937528946 | -0.6051872261780649 | 0.0005603863853470559 | 43251 |

## Interpretation

- `c0` is the irreducible variance (measurement noise).

- `c1` scales the binomial sampling term; it should be stable across landmark sources if the model is correct.

- Four-enzyme and FracMinHash-750 are density-matched (~6,000 landmarks on a 4.5 Mb genome); their MAE and model constants should be most similar.
