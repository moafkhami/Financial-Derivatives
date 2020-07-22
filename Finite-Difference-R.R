E_P_Implicit = function(isCall, K, Tm, 
                       S0, r, sig, N, M, div, dx)
{
  # Implicit Finite Difference Method: i times, 2*i+1 final nodes
  # Precompute constants ----
  #dt = Tm/N
  nu = r - div - 0.5 * sig^2
  edx = exp(dx)
  pu = -0.5 * dt * ( (sig/dx)^2 + nu/dx )
  pm =  1.0 + dt *   (sig/dx)^2 + r*dt 
  pd = -0.5 * dt * ( (sig/dx)^2 - nu/dx)
  firstRow = 1
  nRows = lastRow = 2*M+1
  firstCol = 1
  middleRow = M+1
  nCols = lastCol = N+1
  
  cp = ifelse(isCall, 1, -1)
  # Intialize asset price, derivative price, primed probabilities  ----
  # pp = pmp = V = S = matrix(0, nrow=nRows, ncol=nCols, dimnames=list(
  #   paste("NumUps=",(nCols-1):-(nCols-1), sep=""),
  #   paste("Time=",round(seq(0, 1, len=nCols),4),sep="")))
  pp = pmp = V = S = matrix(0, nrow=nRows, ncol=nCols)
  S[middleRow, firstCol] = S0
  for (i in 1:(nCols-1)) {
    for(j in (middleRow-i+1):(middleRow+i-1)) {
      S[j-1, i+1] = S[j, i] * exp(dx)
      S[j ,  i+1] = S[j, i] 
      S[j+1, i+1] = S[j, i] * exp(-dx)
    }
  }
  # Intialize option values at maturity ----
  for (j in firstRow:lastRow) {
    V[j, lastCol] = max( 0, cp * (S[j, lastCol]-K))
  }
  # Compute Derivative Boundary Conditions ----
  lambdaL = -1 * (S[lastRow-1, lastCol] - S[lastRow,lastCol])
  lambdaU = 0
  
  # Step backwards through the lattice ----
  for (i in (lastCol-1):firstCol) {
    h = solveImplicitTridiagonal(V, pu, pm, pd, lambdaL, lambdaU, i)
    pmp[,i] = h$pmp  # collect the pm prime probabilities
    pp [,i] = h$pp   # collect the p prime probabilities
    V = h$V
    # Apply Early Exercise condition ----
    # for(j in lastRow:firstRow) {
    #   V[j, i] = max(V[j, i], cp * (S[j, lastCol] - K))
    # }
  }
  # Return the price ----
  list(Type = paste( "European", ifelse(isCall, "Call", "Put")),Price = V[middleRow,firstCol],
       Probs=round(c(pu=pu, pm=pm, pd=pd),middleRow), pmp=round(pmp,4), pp=round(pp,4),
       S=round(S,2), V=round(V,middleRow))
}

E_C_Implicit = function(isCall, K, Tm,
                        S0, r, sig, N, M, div, dx)
{
  # Implicit Finite Difference Method: i times, 2*i+1 final nodes
  # Precompute constants ----
  #dt = Tm/N
  nu = r - div - 0.5 * sig^2
  edx = exp(dx)
  pu = -0.5 * dt * ( (sig/dx)^2 + nu/dx )
  pm =  1.0 + dt *   (sig/dx)^2 + r*dt
  pd = -0.5 * dt * ( (sig/dx)^2 - nu/dx)
  firstRow = 1
  nRows = lastRow = 2*M+1
  firstCol = 1
  middleRow = M+1
  nCols = lastCol = N+1

  cp = ifelse(isCall, 1, -1)

  # Intialize asset price, derivative price, primed probabilities  ----
  # pp = pmp = V = S = matrix(0, nrow=nRows, ncol=nCols, dimnames=list(
  #   paste("NumUps=",(nCols-1):-(nCols-1), sep=""),
  #   paste("Time=",round(seq(0, 1, len=nCols),4),sep="")))
  pp = pmp = V = S = matrix(0, nrow=nRows, ncol=nCols)
  S[middleRow, firstCol] = S0
  for (i in 1:(nCols-1)) {
    for(j in (middleRow-i+1):(middleRow+i-1)) {
      S[j-1, i+1] = S[j, i] * exp(dx)
      S[j ,  i+1] = S[j, i]
      S[j+1, i+1] = S[j, i] * exp(-dx)
    }
  }
  # Intialize option values at maturity ----
  for (j in firstRow:lastRow) {
    V[j, lastCol] = max( 0, cp * (S[j, lastCol]-K))
  }
  # Compute Derivative Boundary Conditions ----
  lambdaU = -1 * (S[lastRow-1, lastCol] - S[lastRow,lastCol])
  lambdaL = 0

  # Step backwards through the lattice ----
  for (i in (lastCol-1):firstCol) {
    h = solveImplicitTridiagonal(V, pu, pm, pd, lambdaL, lambdaU, i)
    V = h$V
  }
  # Return the price ----
  list(Type = paste( "European", ifelse(isCall, "Call", "Put")),Price = V[middleRow,firstCol],
       Probs=round(c(pu=pu, pm=pm, pd=pd),middleRow), pmp=round(pmp,4), pp=round(pp,4),
       S=round(S,2), V=round(V,middleRow))
}


