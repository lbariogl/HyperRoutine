o2-analysis-event-selection -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 |
o2-analysis-lf-hypertriton-reco-task -b --configuration json://configuration.json --aod-writer-keep dangling --aod-writer-ntfmerge 300 --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000 |
o2-analysis-timestamp -b --configuration json://configuration.json --shm-segment-size 25000000000 --aod-memory-rate-limit 25000000000
