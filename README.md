# MESS: Multi-Echo Steady-State Free Precession

This is a Python implementation of the MESS (Multi-Echo SSFP) sequence for MRI. 

MESS sequence allows for the simutaneous acquisition of multiple orders of echoes within a single TR.

The code depends on the pypulseq library for pulse sequence design.

```bash
pip install -r requirements.txt
```

To get started, you can run the demo script [`demo.py`](demo.py) to obtain the .seq files for the MESS sequence and bSSFP sequence in the `seq` folder. 


