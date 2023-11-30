# Voltage, V (ger: U)
- Phase Voltage
- Voltage Phase to Ground/Neutral
- Line-to-Line Voltage
- Rated Voltage: The voltage referenced for operational characteristics
- Rated Voltage: Maximum voltage during normal operation >= Rated Voltage

# Current, I
- Rated Current: Maximum current at which an asset can be operated permanently

# Power, S
- Total power of a three-phase system:
  - \( S = 3 \cdot U \cdot I \)
  - \( S = P + jQ \)
  - \( P \): Active Power (=> Heat)
  - \( Q \): Reactive Power (=> Build-Up/Breakdown of Field)

# Impedance, Z
- Impedance: General term for admittance and impedance
- Impedance: Complex resistance \( Z = U/I = R + jX \)
- Resistance: Real part of impedance, \( R \)
  - Ideally constant, neither voltage-dependent nor frequency-dependent
- Reactance: Imaginary part of impedance, \( X \)
  - Ideally only frequency-dependent
  - Coil: \( XL = 2 \pi f L \), \( Z = j \cdot XL \) (where \( f \): Frequency, \( L \): Inductance)
  - Capacitor: \( XC = -1/(2 \pi f C) \), \( Z = -j \cdot XC \) (where \( C \): Capacitance)

# Admittance, Y
- Admittance \( Y \) is the complex conductance (reciprocal of \( Z \), \( Y = 1/Z = 1/(R + jX) \))
  - \( Y = G + jB \)
  
# Frequency
- Amplitude is time-dependent, assumed to be sinusoidal
  - \( u = u_{\text{max}} \cdot \sin(2\pi t + \phi) = u_{\text{max}} \cdot \sin(\omega t + \phi) \)
  - \( \omega \): Angular frequency (2\pi)
  - Note: Angles are in radians!
  - In the transmission network, the frequency is nearly constant at 50Hz
