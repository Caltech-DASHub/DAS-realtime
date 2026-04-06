# DAS-realtime

A Python package to process distributed acousting sensing real-time data streams for earthquake monitoring and early warning. This package can also stream selected channels using a [PyEarthworm](https://github.com/Boritech-Solutions/PyEarthworm) paradigm.

## Installation

In order to install DAS-realtime, you need to have first create a conda environment and install all the required packages. Run the following commands after cloning the repo.

```
conda env create -f environment.yml

git submodule update --init --recursive external/EQNet

git submodule update --init --recursive external/stomp_client

cd external/EQNet/

pip install -r requirements.txt

pip install obspy fastapi

cd ../..

If you want stomp_client to track the latest commit on its configured branch, run:

git submodule sync

git submodule update --remote external/stomp_client

```

In addition, you will need to have PyEarthworm installed as well as Earthworm. Follow the installation guide within the [PyEarthworm](https://github.com/Boritech-Solutions/PyEarthworm) repository.

## PGA-to-FinDer Example

The PGA-to-FinDer path is optional and stays disabled unless `--PGA2FinDer 1` is provided.

Real example command:

```bash
python3 ~/projects/realtime/DAS-realtime/python/StreamProcessor.py \
  --host 162.252.88.51 \
  --port 4900 \
  -strTp OptaSense \
  -wrkint 1.0 \
  --ringbuffer 1.0 \
  --filelength 10.0 \
  -strnRt 1 \
  --xmlmeta ~/projects/realtime/data/input_gitlab/DAS_RidgecrestSouth100km.xml \
  -pch ~/projects/realtime/data/input_gitlab/DAS_RidgecrestSouth100km5000ChPicking_with_PGA2PSRRatio.csv \
  --PGA2FinDer 1 \
  --finderConfig ~/projects/realtime/DAS-realtime/external/stomp_client/stomp_client.cfg \
  --nChPGASmooth 250 \
  --PGAOutput ~/projects/realtime/experiments/output \
  --dataDelayDiag 1
```

Notes for this example:

- The StationXML file defines the sparse channel subset exported to FinDer. Its channel index is the reduced-subset position, not the raw DAS channel id.
- The `-pch` CSV provides the reduced DAS channel subset and the `PGA/PSR-Ratio` conversion factors used to compute pseudo-PGA.
- XML channel index `i` maps directly to row position `i` in the `-pch` CSV subset used for picking and PGA conversion.
- PGA is computed across the full `-pch` channel list, then the XML-selected reduced-subset positions are exported to FinDer.
- `--nChPGASmooth 250` means 125 neighboring channels on the left, 125 on the right, plus the center channel.
- `--PGAOutput` writes a CSV of exported PGA rows and the exact outbound FinDer text payload for validation.
- In this example, `--host/--port` point to the DAS stream, while `--finderConfig` points to the STOMP/FinDer broker configuration; the bundled example config currently targets `localhost`.

# Citation

If you are using this software for your research, please, cite the associated publication:

Biondi, E., Tepp, G., Yu, E., Saunders, J. K., Yartsev, V., Black, M., Watkins, M., Bhaskaran, A., Bhadha, R., Zhan, Z., & Husker, A. L. (2025). Real-time processing of distributed acoustic sensing data for earthquake monitoring operations. Manuscript under review at Seismological Research Letters.