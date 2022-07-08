"""
Placeholder module for real catalog functionality
"""
import os
from glob import glob

import cftime
import xarray as xr

USER = os.environ['USER']

component_system = dict(
    pop='ocn',
    cam='atm',
    cice='ice',
)


def _preprocess_pop(ds):
    tb_var = ds.time.attrs["bounds"]
    time_units = ds.time.units
    calendar = ds.time.calendar

    ds['time'] = cftime.num2date(
        ds[tb_var].mean('d2'),
        units=time_units,
        calendar=calendar,
    )
    ds.time.encoding.update(dict(
        calendar=calendar,
        units=time_units,
    ))
    return ds.set_coords(["KMT", "TAREA"])
    

def _get_assets(case, component, stream, freq, variable_id,years,archive_root,):
    """
    list file assets assuming short term archiver has run
    """
    
    glob_expression = f'{archive_root}/{freq}/{case}.{component}.{stream}.{variable_id}.{years}.nc'
    assets = sorted(glob(glob_expression))
    
    assert assets, f'no files found.\nsearched using: {glob_expression}'

    return assets


def to_dataset_dict(case, component, stream, freq, variable_id, years, archive_root=None,  cdf_kwargs={}):
    
    if archive_root is None:
        archive_root = f'/glade/campaign/cesm/development/bgcwg/projects/hi-res_JRA/cases/g.e22.G1850ECO_JRA_HR.TL319_t13.004/output/ocn/proc/tseries/'
        
    assert os.path.exists(archive_root)
    
    if isinstance(case, str):
        case = [case]
        
    if isinstance(component, str):
        component = [component]
    
    if isinstance(stream, str):
        stream = [stream]
        
    if isinstance(freq, str):
        freq = [freq]        

    if isinstance(variable_id, str):
        variable_id = [variable_id]        
    
    if isinstance(years, str):
        years = [years]        

    component_cdf_kwargs = dict(
        pop=dict(
            coords="minimal",
            combine="by_coords",
            compat="override",
            preprocess=_preprocess_pop,
            decode_times=False,
        ),
    )
                        
    dsets = {}
    for case_i, component_i, stream_i, freq_i, variable_id_i, years_i in zip(case, component, stream, freq, variable_id, years):
        assets = _get_assets(case_i, component_i, stream_i, freq_i, variable_id_i, years_i, archive_root)
        key = f'{case_i}.{component_i}.{stream_i}.{variable_id_i}.{years_i}'
        cdf_kwargs = component_cdf_kwargs[component_i]
        
        dsets[key] = xr.open_mfdataset(assets, **cdf_kwargs)
                
    return dsets