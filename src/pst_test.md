
Ein symmetrischer Phasenschieber (auch als symmetrischer Phasenshifttransformator oder Phase Shifting Transformer (PST) bekannt) wird verwendet, um die Phasenlage und somit den Leistungsfluss in einem Netz zu steuern. Die Y-Parameter eines symmetrischen Phasenschiebers können berechnet werden, indem die Admittanzmatrix des Geräts bestimmt wird.

### Schritte zur Berechnung der Y-Parameter eines symmetrischen Phasenschiebers:

1. **Definition der Basisparameter:**
   - \( Z_{series} \): Serienimpedanz (Widerstand und Reaktanz) des Phasenschiebers.
   - \( Z_{shunt} \): Shuntimpedanz (Widerstand und Reaktanz) des Phasenschiebers.
   - \( \theta \): Phasenverschiebungswinkel.

2. **Bestimmung der Admittanzen:**
   - Serienadmittanz \( Y_{series} = \frac{1}{Z_{series}} \).
   - Shuntadmittanz \( Y_{shunt} = \frac{1}{Z_{shunt}} \).

3. **Admittanzmatrix des Phasenschiebers:**
   Die Y-Parameter werden durch die Kombination der Serien- und Shuntadmittanzen sowie des Phasenverschiebungswinkels berechnet.

### Mathematische Darstellung:

Für einen symmetrischen Phasenschieber, der eine Phasenverschiebung \(\theta\) erzeugt, können die Y-Parameter wie folgt berechnet werden:

- \( Y_{11} = Y_{series} + \frac{1}{2} Y_{shunt} \)
- \( Y_{12} = -Y_{series} \cdot e^{-j\theta} \)
- \( Y_{21} = -Y_{series} \cdot e^{j\theta} \)
- \( Y_{22} = Y_{series} + \frac{1}{2} Y_{shunt} \)

### Beispiel:

Angenommen, wir haben einen symmetrischen Phasenschieber mit den folgenden Parametern:
- Serienimpedanz \( Z_{series} = 0.01 + j0.05 \) Ω.
- Shuntimpedanz \( Z_{shunt} = 0.02 + j0.04 \) Ω.
- Phasenverschiebungswinkel \(\theta = 30^\circ\).

Rechnen wir die Admittanzen und Y-Parameter aus:

```julia
using LinearAlgebra

# Gegebene Parameter
Z_series = Complex(0.01, 0.05)
Z_shunt = Complex(0.02, 0.04)
theta = 30 * (π / 180)  # Umwandlung in Radiant

# Admittanzen berechnen
Y_series = inv(Z_series)
Y_shunt = inv(Z_shunt)

# Y-Parameter berechnen
Y_11 = Y_series + 0.5 * Y_shunt
Y_12 = -Y_series * exp(-im * theta)
Y_21 = -Y_series * exp(im * theta)
Y_22 = Y_series + 0.5 * Y_shunt

# Ergebnisse anzeigen
(Y_11, Y_12, Y_21, Y_22)
```

### Berechnung:

1. **Berechnung der Admittanzen:**

\[
Y_{series} = \frac{1}{0.01 + j0.05} = \frac{1}{0.01 + j0.05} = \frac{0.01 - j0.05}{(0.01)^2 + (0.05)^2} = \frac{0.01 - j0.05}{0.0026} = 3.846 - j19.231 \, \text{S}
\]

\[
Y_{shunt} = \frac{1}{0.02 + j0.04} = \frac{1}{0.02 + j0.04} = \frac{0.02 - j0.04}{(0.02)^2 + (0.04)^2} = \frac{0.02 - j0.04}{0.002} = 10 - j20 \, \text{S}
\]

2. **Berechnung der Y-Parameter:**

\[
Y_{11} = Y_{series} + 0.5 \cdot Y_{shunt} = (3.846 - j19.231) + 0.5 \cdot (10 - j20) = 3.846 - j19.231 + 5 - j10 = 8.846 - j29.231 \, \text{S}
\]

\[
Y_{12} = -Y_{series} \cdot e^{-j\theta} = -(3.846 - j19.231) \cdot e^{-j\frac{\pi}{6}} = -(3.846 - j19.231) \cdot (\cos(-\frac{\pi}{6}) + j\sin(-\frac{\pi}{6})) = -(3.846 - j19.231) \cdot (\frac{\sqrt{3}}{2} - j\frac{1}{2})
\]

\[
Y_{21} = -Y_{series} \cdot e^{j\theta} = -(3.846 - j19.231) \cdot e^{j\frac{\pi}{6}} = -(3.846 - j19.231) \cdot (\cos(\frac{\pi}{6}) + j\sin(\frac{\pi}{6})) = -(3.846 - j19.231) \cdot (\frac{\sqrt{3}}{2} + j\frac{1}{2})
\]

\[
Y_{22} = Y_{11} = 8.846 - j29.231 \, \text{S}
\]

Mit diesem Ansatz können die Y-Parameter eines symmetrischen Phasenschiebers berechnet werden. Der Julia-Code zeigt die Implementierung und Berechnung der entsprechenden Parameter.



