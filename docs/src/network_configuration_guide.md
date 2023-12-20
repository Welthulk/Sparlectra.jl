# Network Configuration Guide

This guide provides detailed explanations for each entry in the provided JSON file describing a power system network.

## Net

- **Value**: String
- **Description**: Unique identifier for the network, e.g. "MyLittleNet"

## BaseMVA

- **Value**: Float
- **Description**: Base MVA for the network data, e.g. 100.0

## Buses

- **Description**: List of buses in the network
  - `bus`: Integer, Bus ID, e.g. 1
  - `name`: String, Bus name, e.g. "B1"
  - `id`: String, Bus identifier, e.g. B1_12345, optional
  - `vn_kv`: Float, Nominal voltage (kV), e.g. 22.0
  - `type`: Integer, Bus type (1: PQ Bus, 2: PV Bus, 3: Slack Bus)
>**Note**:
The order of the bus numbers is arbitrary, but each bus number and ID must be unique. If an ID is not explicitly specified, it will be automatically generated. Additionally, it is required to define only one bus as the slack bus in the network. In a network with three-winding transformers, auxiliary buses are automatically generated.
## Shunts

- **Description**: List of shunts in the network
  - `bus`: Integer, Bus ID, e.g. 1, connection to bus
  - `name`: String, Shunt name, e.g. "Sh1"
  - `id`: String, Shunt identifier, e.g. "Sh1_1234", optional
  - `p_mw`: Float, Active power (MW), e.g. 0.0
  - `q_mvar`: Float, Reactive power (MVAR), e.g. 100.0
  - `vn_kv`: Float, Nominal voltage (kV), e.g. 22.0
  - `step`: Integer, Shunt step, e.g. 1
  - `in_service`: Integer, Shunt in service (1: Yes, 0: No)

## Loads

- **Description**: List of loads in the network
  - `bus`: Integer, Bus ID, e.g. 1, connection to bus
  - `name`: String, Load name, e.g. "L1"
  - `id`: String, Load identifier, e.g. "L1_1234", optional
  - `p_mw`: Float, Active power (MW), e.g. 0.0
  - `q_mvar`: Float, Reactive power (MVAR), e.g. 100.0
  - `in_service`: String, Load in service (1: Yes, 0: No)

## Generators

### StaticGenerators

- **Description**: List of static generators in the network
  - `bus`: String, Bus ID, e.g. 1, connection to bus
  - `name`: String, Generator name, e.g. "G1"
  - `id`: String, Generator identifier, e.g. "G1_1234", optional
  - `p_mw`: Float, Active power (MW), e.g. 10.0
  - `q_mvar`: Float, Reactive power (MVAR), e.g. 20.0
  - `in_service`: String, Generator in service (1: Yes, 0: No)
>**Note**:
In the context of this network configuration, static generators are treated as PQ buses.

### VoltageControlledGenerators

- **Description**: List of voltage-controlled generators in the network, each represented by a dictionary with the following properties:
  - `bus`: String, Bus ID, e.g. 1, connection to bus
  - `name`: String, Generator name, e.g. "G1"
  - `id`: String, Generator identifier, e.g. "G1_1234", optional
  - `p_mw`: Float, Active power (MW), e.g. 10.0
  - `max_q_mvar`: Float, Maximum reactive power (MVAR)
  - `min_q_mvar`: Float, Minimum reactive power (MVAR)
  - `vm_pu`: Float, Voltage magnitude (pu), e.g. 1.05
  - `in_service`: Integer, Generator in service (1: Yes, 0: No)
>**Note**:
`max_q_mvar`: not used in the current version of the code
`min_q_mvar`: not used in the current version of the code
In the context of this network configuration, static generators are treated as PV buses.

### ExternalGrids
- **Description**: List of external grids in the network
  - `bus`: String, Bus ID, e.g. 1, connection to bus
  - `name`: String, Generator name, e.g. "Ext1"
  - `id`: String, Generator identifier, e.g. "Ext1_1234", optional  
  - `vm_pu`: Float, Voltage magnitude (pu), e.g. 1.05
  - `va_degree`: Float, Voltage angle (degree), e.g. 0.0
  - `in_service`: Integer, Generator in service (1: Yes, 0: No)
>**Note**:
In the context of this network configuration, external grids are designated as slack buses, allowing for phase angle specification.



## LineModelling

- **Description**: List of line models in the network
  - `model`: String, Line model, e.g. "STD1"
  - `id`: String, Line model identifier, e.g. "STD1_1234", optional
  - `r_ohm_per_km`: Float, Resistance per kilometer (ohm/km), e.g. 0.2
  - `x_ohm_per_km`: Float, Reactance per kilometer (ohm/km), e.g. 0.39
  - `c_nf_per_km`: Float, Capacitance per kilometer (nF/km),  e.g. 9.55
  - `g_us_per_km`: Float, Conductance per kilometer (uS/km), e.g. 0.0
  - `max_i_ka`: Float, Maximum current (kA), e.g. 1000.0
>**Note**:
`max_i_ka`: not used in the current version of the code
### ACLines using LineModelling

- **Description**: List of AC lines in the network
  - `from_bus`: Integer, Starting bus ID, e.g. 1
  - `to_bus`: Integer, Ending bus ID, e.g. 2
  - `name`: String, Line name, e.g. "L1_2"
  - `id`: String, Line identifier, e.g. "L1_2_1234", optional
  - `length_km`: Float, Line length (km), e.g. 1.0
  - `line_model`: Integer, Line model, e.g. "STD1"
  - `parallel`: Integer, Number of parallel lines, e.g. 1, optional
  - `in_service`: Integer, Line in service (1: Yes, 0: No)
>**Note**:
`parallel`: ENTER TEXT HERE
### ACLines
 - **Description**: List of AC lines in the network
   - `from_bus`: Integer, Starting bus ID, e.g. 1
   - `to_bus`: Integer, Ending bus ID, e.g. 2
   - `name`: String, Line name, e.g. "L1_2"
   - `id`: String, Line identifier, e.g. "L1_2_1234", optional
   - `length_km`: Float, Line length (km), e.g. 1.0
   - `r_ohm_per_km`: Float, Resistance per kilometer (ohm/km), e.g. 0.2
   - `x_ohm_per_km`: Float, Reactance per kilometer (ohm/km), e.g. 0.39
   - `c_nf_per_km`: Float, Capacitance per kilometer (nF/km),  e.g. 9.55
   - `g_us_per_km`: Float, Conductance per kilometer (uS/km), e.g. 0.0
   - `max_i_ka`: Float, Maximum current (kA), e.g. 1000.0
   - `in_service`: Integer, Line in service (1: Yes, 0: No)

## TransformerModelling
### VK-Characteristics

- **Description**: Dictionary of voltage characteristics of tap position for transformers, each represented by a key-value pair with the model name as the key and a dictionary as the value containing `tap` and `vk` arrays. 
  - `tap`: Array of tap positions, e.g. [-9, 0, 9]
  - `vk`: Array of voltage magnitudes, e.g. [19, 12, 19]
>**Note**:
A spline function is created using the provided pair of values, `tap`, and `vk`. This spline function serves to determine the voltage dependency of the equivalent circuit elements based on the tap position.
A voltage characteristic model is optional. 
### TapChangerModelling

- **Description**: List of tap changer models in the network
  - `model`: String, Tap changer model, e.g. "VW3"
  - `id`: String, Tap changer model identifier, e.g. "VW3_1234", optional
  - `tap_side`: String, Tap side, e.g., "HV", allowed values: "HV", "MV", "LV"
  - `tap_min`: Integer, Minimum tap position, e.g. -9
  - `tap_max`: Integer, Maximum tap position, e.g. 9
  - `tap_neutral`: Integer, Neutral tap position, e.g. 0
  - `tap_step_percent`: Float, Tap step percentage, e.g. 1.25
  - `neutralU_ratio`: Float, Neutral U ratio, e.g. 0.8
  - `tap_step_degree`: Float, Tap step degree, e.g. 10.0
  - `shift_degree`: Float, Shift degree, e.g. 0.0

### TwoWindingTransformers

- **Description**: List of two-winding transformers in the network
  - `hv_bus`: Integer, High-voltage bus ID, e.g. 1
  - `lv_bus`: Integer, Low-voltage bus ID, e.g. 2
  - `name`: String, Transformer name, e.g. "WT2_12"
  - `id`: String, Transformer identifier, e.g. "WT2_12_1234", optional
  - `vn_hv_kv`: Float, High-voltage nominal voltage (kV), e.g. 220.0
  - `vn_lv_kv`: Float, Low-voltage nominal voltage (kV), e.g. 110.0
  - `vk_percent`: Float, Short-circuit voltage (%), e.g. 10.0
  - `vkr_percent`: Float, Short-circuit voltage (%), e.g. 0.41
  - `pfe_kw`: Float, Iron losses (kW), e.g. 12.0
  - `i0_percent`: Float, No-load losses (%), e.g. 0.007
  - `tap_model`: Float, Tap changer model, e.g. "VW3"
  - `tap_pos`: Integer, Tap position, e.g. 3
  - `vk_dependence`: String, Voltage characteristic model, e.g. "VK1", optional
  - `in_service`: Transformer in service (1: Yes, 0: No)

### ThreeWindingTransformers

- **Description**: List of three-winding transformers in the network
  - `hv_bus`: Integer, High-voltage bus ID, e.g. 1
  - `mv_bus`: Integer, Medium-voltage bus ID, e.g. 2
  - `lv_bus`: Integer, Low-voltage bus ID, e.g. 3
  - `name`: String, Transformer name, e.g. "WT3_123"
  - `id`: String, Transformer identifier, e.g. "WT3_123_1234", optional
  - `vn_hv_kv`: Float, High-voltage nominal voltage (kV), e.g. 220.0
  - `vn_mv_kv`: Float, Medium-voltage nominal voltage (kV), e.g. 110.0
  - `vn_lv_kv`: Float, Low-voltage nominal voltage (kV), e.g. 30.0
  - `sn_hv_mva`: Float, High-voltage nominal power (MVA), e.g. 250.0
  - `sn_mv_mva`: Float, Medium-voltage nominal power (MVA), e.g. 150.0
  - `sn_lv_mva`: Float, Low-voltage nominal power (MVA), e.g. 50.0
  - `vk_hv_percent`: Float, High-voltage short-circuit voltage (%), e.g. 10.0
  - `vk_mv_percent`: Float, Medium-voltage short-circuit voltage (%), e.g. 11.0
  - `vk_lv_percent`: Float, Low-voltage short-circuit voltage (%), e.g. 12.0
  - `vkr_hv_percent`: Float, High-voltage short-circuit voltage (%), e.g. 0.3
  - `vkr_mv_percent`: Float, Medium-voltage short-circuit voltage (%), e.g. 0.31
  - `vkr_lv_percent`: Float, Low-voltage short-circuit voltage (%), e.g. 0.32
  - `pfe_kw`: Float, Iron losses (kW), e.g. 30.0
  - `i0_percent`: Float, No-load losses (%), e.g. 0.1
  - `shift_mv_degree`: Float, Phase shift (degree), e.g. 0.0
  - `shift_lv_degree`: Float, Phase shift (degree), e.g. 0.0
  - `tap_model`: Integer, Tap changer model, e.g. 2
  - `tap_pos`: Integer, Tap position, e.g. 9
  - `vk_dependence`: String, Voltage characteristic model, e.g. "VK1", optional
  - `in_service`: Integer, Transformer in service (1: Yes, 0: No)


  # Example Network
  
  ## bsp6, 4-Bus Test Network, 110kV, mBase:100.0 MVA


              4                  1
       ExtG1->|-------L4-1-------|<-SG1
              |                  
              |
              |                  2            3 (22kV)                   
              |-------L4-2-------|------T1----|->LD1
                            VG1->|            x Sh1
```json
{    
    "Net": "bsp6",
    "BaseMVA": 100.0,
    "Buses": [
      {"bus": 1, "name": "Bus_1", "id": "", "vn_kv": 110, "type": 1},
      {"bus": 2, "name": "Bus_2", "id": "", "vn_kv": 110, "type": 3},
      {"bus": 3, "name": "Bus_3", "id": "", "vn_kv": 22,  "type": 1},
      {"bus": 4, "name": "Bus_4", "id": "", "vn_kv": 110, "type": 2}
      ],
    "Shunts": [ 
      { "bus": 3, "name": "SH1", "id": "", "p_mw": 0.0, "q_mvar": 3.8, "vn_kv":22, "step": 1, "in_service": 1 }            
    ],
    "Loads": [
      { "bus": 3, "name": "LD1", "id": "", "p_mw": 2.0, "q_mvar": 4.0, "in_service": 1 }
    ],
    "StaticGenerators": [
      { "bus": 1, "name": "SG1", "id": "", "p_mw": 2.0, "q_mvar": -0.5, "in_service": 1 }
    ],
    "VoltageControlledGenerators": [
      { "bus": 2, "name": "VG1", "id": "", "p_mw": 6.0, "max_q_mvar": 15.0,  "min_q_mvar": -15.0, "vm_pu": 1.03, "in_service": 1 }
    ],
    "ExternalGrids": [
      { "bus": 4,"name": "ExG1", "id": "", "vm_pu": 1.02, "va_degree": 0.0, "in_service": 1  }
    ],
    "LineModelling":[
      {"model":"STD1", "id": "", "r_ohm_per_km": 0.2, "x_ohm_per_km": 0.39, "c_nf_per_km": 9.55, "g_us_per_km": 0.0, "max_i_ka": 0.0}
    ],
    
    "ACLines": [
      { "from_bus": 4,"to_bus": 2, "name": "L4_2", "id": "", "length_km": 25.0, "line_model": "STD1", "in_service": 1 },
      { "from_bus": 4,"to_bus": 1, "name": "L4_1","id": "", "length_km": 10.0, "r_ohm_per_km": 0.05, "x_ohm_per_km": 0.39, "c_nf_per_km": 9.55, "g_us_per_km": 0.0, "max_i_ka": 0.0, "in_service": 1}
    ],
    "VK-Characteristics":{
       "model1":{
         "tap":[1,    4,    9],
         "vk":[7.0,  10.0,  12.0]},
       "model2":{
          "tap":[-9, 0, 9],
          "vk":[19, 12, 19]}
    },    
    "TapChangerModelling": [
      { "model": 1, "id": "", "tap_side": "HV", "tap_min": 1, "tap_max": 9,  "tap_neutral": 4, "tap_step_percent": 2.5, "neutralU_ratio": 1.0, "tap_step_degree": 0.0, "shift_degree": 0.0 }
    ],
    "TwoWindingTransformers": [
      {"hv_bus": 2,"lv_bus": 3,"name": "T1", "id": "", "sn_mva": 25.0, "vn_hv_kv": 110.0, "vn_lv_kv": 22.0, "vk_percent": 12.0, "vkr_percent": 0.41, "pfe_kw": 14.0, "i0_percent": 0.07, "tap_model":1, "tap_pos": 4, "vk_dependence": "model1", "in_service": 1}
      ] 
}
```
# Parsing Function
The function `createNetFromFile` to parse the JSON file is implemented in the include-file [`readnetfromfile.jl`](../../src/readnetfromfile.jl). The function takes as input the filepath to the JSON file and returns a network object. The function is implemented as follows:
```js
function createNetFromFile(filename, base_MVA::Float64 = 0.0, log::Bool = false)::ResDataTypes.Net
```
>**Note**:
The `base_MVA` parameter is optional and defaults to 0.0. If the `base_MVA` parameter is not provided, the base MVA is read from the JSON file. If the `base_MVA` parameter is provided, the base MVA from the JSON file is ignored.