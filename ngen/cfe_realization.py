#!/usr/bin/env python3


import json
import pandas
import typing
import pathlib
from collections import OrderedDict

class GlobalRealization:
    def __init__(self, global_params=None, time=None, routing=None, catchment_realizations=None):
        self.global_params = global_params
        self.time = time
        self.routing = routing
        self.catchment_realizations = catchment_realizations
    def toJSON(self):
        
        # construct the object that we want to return
        dat = {}
        if self.global_params is not None:
            dat['global'] = self.global_params
        if self.catchment_realizations is not None:
            dat.update(json.loads(self.catchment_realizations.toJSON()))
        if self.time is not None:
            dat['time'] = self.time
        if self.routing is not None:
            dat['routing'] = self.routing
        
        return json.dumps(dat, sort_keys=False, indent=4)

class CatchmentRealizations:
    def __init__(self):
        self.catchments = {}
    def add_realization(self, name, realization):
        self.catchments[name] = realization
    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, 
            sort_keys=False, indent=4)

class Realization:
    def __init__(self, formulation, forcing={}):
        self.formulations = [formulation]
        self.forcing = forcing
    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, 
            sort_keys=False, indent=4)

class Formulation:
    def __init__(self, name, params={}, modules=[]):
        self.name = name
        self.params = params
        self.params['modules'] = modules
    def add_module(self, module):
        self.modules.append(module)
    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, 
            sort_keys=False, indent=4)
        
class Module:
    def __init__(self, name, params={}):
        self.name = name
        self.params=params
    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, 
            sort_keys=False, indent=4)
    
    
def parse_cfe_parameters(cfe_noahowp_attributes: pandas.DataFrame) -> typing.Dict[str, dict]:
    """
    Parses parameters from NOAHOWP_CFE DataFrame
    
    Parameters
    ----------
    cfe_noahowp_attributes: pandas.DataFrame
        Dataframe of NoahOWP CFE parameters
    
    Returns
    -------
    Dict[str, dict]: parsed CFE parameters
    
    """
    catchment_configs = {}
    for idx, row in cfe_noahowp_attributes.iterrows():
        d = OrderedDict()

        # static parameters
        d['forcing_file']='BMI'
        d['surface_partitioning_scheme']='Schaake'

        # ----------------    
        # State Parameters
        # ----------------

        # soil depth
        d['soil_params.depth']='2.0[m]'

        # many of these values are taken from the 2m depth in hydrofabrics cfe_noahowp_attributes
        d['soil_params.b']=f'{row["bexp_soil_layers_stag=2"]}[]' # 	beta exponent on Clapp-Hornberger (1978) soil water relations

        # saturated hydraulic conductivity
        d['soil_params.satdk']=f'{row["dksat_soil_layers_stag=2"]}[m s-1]' 

        # saturated capillary head
        d['soil_params.satpsi']=f'{row["psisat_soil_layers_stag=2"]}[m]'

        # this factor (0-1) modifies the gradient of the hydraulic head at the soil bottom. 0=no-flow.
        d['soil_params.slop']=f'{row["slope"]}[m/m]'

        # saturated soil moisture content
        d['soil_params.smcmax']=f'{row["smcmax_soil_layers_stag=2"]}[m/m]'

        # wilting point soil moisture content
        d['soil_params.wltsmc']=f'{row["smcwlt_soil_layers_stag=2"]}[m/m]'

        # ---------------------    
        # Adjustable Parameters
        # ---------------------

        # optional; defaults to 1.0
        d['soil_params.expon']=f'{row["gw_Expon"]}[]' if row["gw_Expon"] is not None else '1.0[]'

        # optional; defaults to 1.0 
        # not sure if this is the correct key
        d['soil_params.expon_secondary']=f'{row["gw_Coeff"]}[]' if row["gw_Coeff"] is not None else '1.0[]'

        # maximum storage in the conceptual reservoir
        d['max_gw_storage']=f'{row["gw_Zmax"]}[m]' if row["gw_Zmax"] is not None else '0.011[m]'

         # primary outlet coefficient
        d['Cgw']='0.0018[m h-1]'

        # exponent parameter (1.0 for linear reservoir)
        d['expon']='6.0[]' 

        # initial condition for groundwater reservoir - it is the ground water as a 
        # decimal fraction of the maximum groundwater storage (max_gw_storage) for the initial timestep
        d['gw_storage']='0.05[m/m]'

        # field capacity
        d['alpha_fc']='0.33' 

        # initial condition for soil reservoir - it is the water in the soil as a 
        # decimal fraction of maximum soil water storage (smcmax * depth) for the initial timestep
        d['soil_storage']='0.05[m/m]' 

        # number of Nash lf reservoirs (optional, defaults to 2, ignored if storage values present)
        d['K_nash']='0.03[]' 

        # Nash Config param - primary reservoir
        d['K_lf']='0.01[]' 

        # Nash Config param - secondary reservoir
        d['nash_storage']='0.0,0.0' 

        # Giuh ordinates in dt time steps
        d['giuh_ordinates']='1.00,0.00'

        # ---------------------    
        # Time Info
        # ---------------------


        # set to 1 if forcing_file=BMI
        d['num_timesteps']='1' 

        # ---------------------    
        # Options
        # ---------------------

        # prints various debug and bmi info
        d['verbosity']='1' 

        d['DEBUG']='0'

        # Parameter in the surface runoff parameterization 
        # (https://mikejohnson51.github.io/hyAggregate/#Routing_Attributes)
        d['refkdt']=f'{row["refkdt"]}'

        catchment_configs[row.divide_id] = d
    
    return catchment_configs

def make_catchment_configs(base_dir: pathlib.Path, catchment_configs: pandas.DataFrame) -> None:
    """
    Creates CFE configuration files for each catchment in "catchment_configs"
    
    Parameters
    ----------
    base_dir: pathlib.Path
        base directory to save the output configuration files
    catchment_configs: pandas.DataFrame
        dataframe containing NoahOWP CFE parameters
        
    Returns
    -------
    None
    
    """
    
    for name, conf in catchment_configs.items():
        with open(f'{base_dir}/{name}_config.ini', 'w') as f:
            for k, v in conf.items():
                f.write(f'{k}={v}\n')

def create_cfe_realization(base_dir: pathlib.Path,
                           cfe_noahowp_csv: pathlib.Path,
                           time={},
                           config_path=pathlib.Path('.'),
                           forcing_path=pathlib.Path('.'),
                           troute_path=None,
                           binary_path=pathlib.Path('/dmod/shared_libs'),
                           ):
    
    catchment_configs = parse_cfe_parameters(pandas.read_csv(cfe_noahowp_csv))
    make_catchment_configs(base_dir, catchment_configs)
    
    catchment_realizations = CatchmentRealizations()

    for key, val in catchment_configs.items():
        config_name = f'{key}_config.ini'
        config_path_ini = f'{config_path}/{config_name}'
        forcing_file_path = f'{forcing_path}/{key}.csv'

        # CFE
        module_params = {"name": "bmi_c",
                         "model_type_name": "CFE",
                         "main_output_variable": "Q_OUT",
                         "init_config": f"{config_path_ini}",
                         "allow_exceed_end_time": True,
                         "fixed_time_step": False,
                         "uses_forcing_file": True,
                         "forcing_file": forcing_file_path,
                         "variables_names_map": {
                             "atmosphere_water__liquid_equivalent_precipitation_rate": "precip_rate",
                             "water_potential_evaporation_flux": "EVAPOTRANS",
                             "ice_fraction_schaake": "sloth_ice_fraction_schaake",
                             "ice_fraction_xinan": "sloth_ice_fraction_xinan",
                             "soil_moisture_profile": "sloth_smp"
                         },
                        "library_file": f"{binary_path}/libcfebmi.so.1.0.0",
                        "registration_function": "register_bmi_cfe"}
        m1 = Module('bmi_c', params=module_params)

        # SLOTH
        module_params = {"name": "bmi_c++",
                         "model_type_name": "SLOTH",
                         "main_output_variable": "z",
                         "init_config": "/dev/null",
                         "allow_exceed_end_time": True,
                         "fixed_time_step": False,
                         "uses_forcing_file": False,
                         "model_params": {
                             "sloth_ice_fraction_schaake(1,double,m,node)": "0.0",
                             "sloth_ice_fraction_xinan(1,double,1,node)": "0.0",
                             "sloth_smp(1,double,1,node)": "0.0",
                             "EVAPOTRANS": "0.0"
                         },
                        "library_file": f"{binary_path}/libslothmodel.so",
                        "registration_function": "none",}
        m2 = Module('bmi_c++', params=module_params)


        form_params = {"name": "bmi_multi",
                       "model_type_name": "NoahOWP_CFE",
                       "main_output_variable": "Q_OUT",
                       "init_config": "",
                       "allow_exceed_end_time": False,
                       "fixed_time_step": False,
                       "uses_forcing_file": False}
        f = Formulation('bmi_multi', params=form_params, modules=[m1, m2])

        realization = Realization(f, forcing={'path': forcing_file_path})

        catchment_realizations.add_realization(key, realization)


    if troute_path is not None:
        # routing = {'t_route_config_file_with_path': '/ngen/ngen/data/config/ngen.yaml'}
        routing = {'t_route_config_file_with_path': f'{troute_path}'}
    else:
        routing = None
    
    realization = GlobalRealization(time=time,
                                    catchment_realizations=catchment_realizations,
                                    routing=routing)
    
    with open(f'{base_dir}/realization.json', 'w') as f:
        f.write(realization.toJSON())
        
        
    ####### ngen-yaml
    import ruamel.yaml
    # load
    yaml = ruamel.yaml.YAML()
    with open('./ngen-routing-template.yaml') as file:
        ngen = yaml.load(file)
    # define wb_id    
    wb_id = base_dir.parts[-2]
    # modify
    ngen['network_topology_parameters']['supernetwork_parameters']['geo_file_path'] = f'{config_path}/{wb_id}_upstream_subset.gpkg'
    ngen['network_topology_parameters']['waterbody_parameters']['level_pool']['level_pool_waterbody_parameter_file_path'] = f'{config_path}/{wb_id}_upstream_subset.gpkg'
    ngen['network_topology_parameters']['waterbody_parameters']['level_pool']['reservoir_parameter_file'] = f'{config_path}/{wb_id}_upstream_subset.gpkg'
    ngen['compute_parameters']['restart_parameters']['start_datetime'] = time['start_time']  #'2022-08-24 13:00:00'
    ngen['compute_parameters']['forcing_parameters']['nts'] = time['nts'] 
    
    # save
    with open(f'{wb_id}/config/ngen.yaml', 'w') as file:
        yaml.dump(ngen, file)
