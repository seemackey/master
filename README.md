# master

Guide to seemackey repository on github

1) "CP" scripts for choice probability/detect probability analysis from
Mackey et al. (2022)

  Preprint  here:
  https://doi.org/10.1101/2022.06.15.496306
  
  J Neurophysiol paper here:
  https://doi.org/10.1152/jn.00439.2022

  Uses: CP_Analysis_ScriptCM,
        CP_Calc,
        CP_ROC

  Shell script is CP_Analysis_ScriptCM.m imports data and runs CP_Calc in a loop
  through all the data. CP_calc runs the CP function, which calculates CP using
  signal detection theoretic ROC analysis (CP_ROC function)


2) "MFdisc" and "samsmodlong" scripts for neuronal discrimination analysis
and population modeling using single-unit responses to amplitude modulated noise

  Preprint: https://doi.org/10.1101/2022.08.05.502987

  Uses: samsmodlong_pooler_looper,
        samsmodlong_pooler_fxn,
        mfSort,
        MTF.m,
        mfROC or mfROC_pop,
        and (!!OPTIONAL!!) trimmer.m if the duration is less than 500 ms, which can trim response to nearest AM trough.
        some reviewers hated trimmer.m. - not currently using it.

  Shell script is samsmodlong_pooler_looper.m, which imports data and runs
  samsmodlong_pooler_fxn - this function uses mfSort.m and MTF.m to pick out
  neurons and certain trials. At that point the responses are extracted and it
  does ROC analysis (mfROC for a single neuron or mfROCpop for populations)
  This is describing the population level analysis, but samsmodlong_noActX.m
  does all of this for a single neuron.
  
  
  3) Drift diffusion model (jupyter notebook) mostly uses code from de Gee et al. (2021, eLife)
  but some was repurposed. Links to de Gee et al's code and papers are in the jupyter notebook "Chase_discrimination_DDM_clean2"
  

