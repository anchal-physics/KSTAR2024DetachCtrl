# Data

This data has been downloaded from KSTAR mds server using the python package [mdsh5](https://pypi.org/project/mdsh5/). To download your own copy or update data with changes in the associated [KSTAR_channels.yml](KSTAR_channels.yml), do following:

```
pip install mdsh5
read_mds -c KSTAR_channels.yml
```