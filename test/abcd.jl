# Funktion zur Berechnung der Admittanzmatrix aus r_pu, x_pu und b_pu basierend auf dem PI-Modell
function calculate_Y_from_PI_model(r_pu, x_pu, b_pu)
    # Series Admittance ys
    ys = inv((r_pu + x_pu*im))
    
    # Shunt Admittance ysh
    ysh = im * b_pu
    
    # Berechne die Admittanzmatrix Y
    Y_from_from = ys + 0.5*ysh
    Y_from_to = -ys
    Y_to_from = -ys
    Y_to_to = ys + 0.5*ysh
    
    Y = [Y_from_from Y_from_to;
         Y_to_from Y_to_to]
        
    return Y
end

# Funktion zur Berechnung der ABCD-Parameter aus r_pu, x_pu und b_pu
function getABCD_from_parameters(r_pu, x_pu, b_pu)
  # Series Admittance ys
    ys = inv(r_pu + x_pu * im)
    
    # Shunt Admittance ysh
    ysh = im * b_pu
    
    # Calculate Y_from_from, Y_from_to, Y_to_from, Y_to_to
    Y_from_from = ys + 0.5 * ysh
    Y_from_to = -ys
    Y_to_from = -ys
    Y_to_to = ys + 0.5 * ysh
    
    # Calculate ABCD parameters
    A = Y_from_from
    B = Y_from_to
    C = Y_to_from
    D = Y_to_to
    
    return A, B, C, D

 end   

# Funktion zur Berechnung der Admittanzmatrix Y aus ABCD-Parametern
function calculate_Y_from_ABCD(A, B, C, D)
    Y_from_from = A
    Y_from_to = B
    Y_to_from = C
    Y_to_to = D
    
    Y = [Y_from_from Y_from_to;
         Y_to_from Y_to_to]
    
    return Y
end


# Beispielwerte für r_pu, x_pu und b_pu
r_pu = 0.1
x_pu = 0.2
b_pu = 0.01

# Berechne ABCD-Parameter aus den gegebenen Werten
A, B, C, D = getABCD_from_parameters(r_pu, x_pu, b_pu)

# Berechne die Admittanzmatrix aus ABCD-Parametern
Y_from_ABCD = calculate_Y_from_ABCD(A, B, C, D)

# Berechne die Admittanzmatrix aus dem PI-Modell zur Überprüfung
Y_from_PI_model = calculate_Y_from_PI_model(r_pu, x_pu, b_pu)

# Vergleiche die beiden Matrizen
println("Admittanzmatrix Y aus PI-Modell:")
println(Y_from_PI_model)
println("Admittanzmatrix Y aus ABCD-Parametern:")
println(Y_from_ABCD)

# Überprüfe, ob die Matrizen gleich sind (mit einer kleinen Toleranz aufgrund von Gleitkommadarstellung)
if isapprox(Y_from_PI_model, Y_from_ABCD; atol=1e-10)
    println("Die Admittanzmatrizen sind identisch.")
else
    println("Fehler: Die Admittanzmatrizen sind nicht identisch.")
end
