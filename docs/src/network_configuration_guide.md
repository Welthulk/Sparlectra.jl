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
The order of the bus numbers is arbitrary, but each bus number and ID must be unique. If an ID is not explicitly specified, it will be automatically generated. Additionally, it is required to define only one bus as the slack bus in the network.
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

# Generators
>**Note**:

## StaticGenerators

- **Description**: List of static generators in the network
  - `bus`: String, Bus ID, e.g. 1, connection to bus
  - `name`: String, Generator name, e.g. "G1"
  - `id`: String, Generator identifier, e.g. "G1_1234", optional
  - `p_mw`: Float, Active power (MW), e.g. 10.0
  - `q_mvar`: Float, Reactive power (MVAR), e.g. 20.0
  - `in_service`: String, Generator in service (1: Yes, 0: No)
>**Note**:
In the context of this network configuration, static generators are treated as PQ buses.

## VoltageControlledGenerators

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

## ExternalGrids
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
## ACLines using LineModelling

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
## ACLines
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

## VK-Characteristics

- **Description**: Dictionary of voltage characteristics models for transformer, each represented by a key-value pair with the model name as the key and a dictionary as the value containing `tap` and `vk` arrays. TODO: ENTER TEXT HERE
  - `tap`: Array of tap positions, e.g. [-9, 0, 9]
  - `vk`: Array of voltage magnitudes, e.g. [19, 12, 19]

## TapChangerModelling

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

## TwoWindingTransformers

- **Description**: List of two-winding transformers in the network (currently empty).

## ThreeWindingTransformers

- **Description**: List of three-winding transformers in the network, each represented by a dictionary with the following properties:
  - `hv_bus`: High-voltage bus ID
  - `mv_bus`: Medium-voltage bus ID
  - `lv_bus`: Low-voltage bus ID
  - `name`: Transformer name
  - `id`: Transformer identifier
  - `vn_hv_kv`: High-voltage nominal voltage (kV)
  - `vn_mv_kv`: Medium-voltage nominal voltage (kV)
  - `vn_lv_kv`: Low-voltage nominal voltage (kV)
  - `sn_hv_mva`: High-voltage nominal power (MVA)
  - `sn_mv_mva`: Medium-voltage nominal power (MVA)
  - `sn_lv_mva`: Low-voltage nominal power (MVA)
  - Other transformer parameters...
