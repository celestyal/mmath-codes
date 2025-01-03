# MMath codes 
## Introduction
This repository contains codes used in my final year Applied Mathematics MMath project.

## How to use the IDL codes
Clone the repository and add the `src` directory to the IDL path. The `tools/` directory contains `update_idl_path` to automate this for the current shell. This should work for most shells. IDL sessions then called will be able to access functions and procedures in `src`.

You can do this with

```console
  git clone https://github.com/celestyal/mmath-codes.git
  cd mmath-codes
  . tools/update_idl_path
```

**NOTE:** By default, `update_idl_path` assumes you are using `gdl`. To change this, manually edit the `IDL_CMD` variable at the start of the script.

## Example codes
Example codes will be stored in `examples/`. Currently the only code example is `kitt_peak.pro`. This script was written to develop an understanding of solar magnetic fields in research, by working with synoptic magnetogram maps collected by the NSO Vacuum Telescope located on Kitt Peak, Arizona.

## Acknowledgment
NSO/Kitt Peak data used here are produced cooperatively by NSF/NOAO, NASA/GSFC, and NOAA/SEL.

The data can be accessed at https://nispdata.nso.edu/ftp/kpvt/synoptic/
