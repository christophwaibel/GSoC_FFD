# FFD-lib
FFD library with parallel loops. required for GH_Wind

Please refer to this publication for citation: [Waibel et al. (2017)](https://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0ahUKEwiJmpLqiN_WAhVHQBoKHUO5DEUQFgg0MAA&url=https://www.conftool.pro/bs2017/index.php/BS2017_Airflow_03_3_2582_Waibel_2017-04-19_03-30_a.pdf?page=downloadPaper%26filename=BS2017_Airflow_03_3_2582_Waibel_2017-04-19_03-30_a.pdf%26form_id=2582%26form_version=final&usg=AOvVaw3kMQmDlXdkID-rco11CyPk)

Modified after: https://github.com/lukasbystricky/GSoC_FFD


To Do:
- [ ] bug: pressure field keeps updating in an open domain
- [ ] add energy transportation (Boussinesque approximation)
- [ ] add pollution transportation (Fick's law)
- [ ] add numeric convergence indicators
- [ ] try different PDE solvers (multigrid, conjugate gradient, TDMA...?) 
<br></br>
- [ ] parallelize on GPU
<br></br>
- [ ] meta models for spatially resolved effective viscosity 
