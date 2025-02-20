# Data

This data has been downloaded from KSTAR mds server using the python package [mdsh5](https://pypi.org/project/mdsh5/). To download your own copy or update data with changes in the associated [KSTAR_channels.yml](KSTAR_channels.yml), do following (replace <KSTARmdsUsername> with your KSTAR MDS username. Ex: guptaa):

```
pip install mdsh5
read_mds -c KSTAR_channels.yml -s <KSTARmdsUsername>@mdsr2.science.kstar.kfe.re.kr:8005
read_mds -c KSTAR_EFIT_Psi.yml -s <KSTARmdsUsername>@mdsr2.science.kstar.kfe.re.kr:8005
```
