This python library was developed specifically to analyze 2-dimensional image Raman scattering data on single crystal and electronic semiconductor materials but can be further used for other type of samples as well. The software was structured so that additional built-in classes and functions can readily be added for other type of analyses. The software utilizes the wdfReader library. The authors would like to thank lovaulonze for creating the wdfReader library (under the open-sourced MIT license)

Note: 

Changes have been made to the wdfreader (https://github.com/alchem0x2A/py-wdf-reader). Possibly need to install renishawWiRE
and cut/paste the below line to line 26 in pyRaman.py file:

from renishawWiRE import WDFReader as wdf
