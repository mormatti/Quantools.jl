# This file contains the functions used to generate the spin matrices. The function
# spinOperator returns the spin operator corresponding to the spin s and the symbol
# corresponding to the operator. For example, spinOperator("1/2", "Sx") returns the
# spin operator S^x for the spin 1/2. The symbol corresponding to the operator can be
# "S+", "S-", "Sx", "iSy", or "Sz".

# We import the package ITensors
using ITensors

"""This function converts a string which represents the spin to an integer which represents
the spin doubled. For example, "1/2" is converted to 1, "1" is converted to 2, and "3/2"
is converted to 3."""
function doubledSpin(s)
  # If s is a string composed of a single character, we convert that single character
  # to an integer n, and we return the double of n.
  if length(s) == 1
    n = parse(Int, s)
    return 2*n
  # If s is a string composed of three characters, we convert the first character to an
  # integer n, and we return n.
  elseif length(s) == 3
    n = parse(Int, s[1])
    return n
  # otherwise, we throw an error.
  else
    error("Invalid spin")
  end
end

"""This function returns the spin operator corresponding to the spin s and the symbol
corresponding to the operator. For example, spinOperator("1/2", "Sx") returns the
spin operator S^x for the spin 1/2. The symbol corresponding to the operator can be
"S+", "S-", "Sx", "iSy", or "Sz"."""
function spinOperator(s, operator)
  # We convert the spin s to an integer which represents the spin doubled.
  n = doubledSpin(s)
  j = n/2
  # We initialize the array which will contain the spin ladder operator.
  result = zeros(n+1, n+1)
  # If direction = "+", we fill the array with the corresponding values.
  if operator == "S+"
    for i in 1:n
      m = i - j - 1
      result[i, i+1] = sqrt(j*(j+1) - m*(m+1))
    end
  # If direction = "-", we fill the array with the corresponding values.
  elseif operator == "S-"
    for i in 2:n+1
      m = i - j - 1
      result[i, i-1] = sqrt(j*(j+1) - m*(m-1))
    end
  # If axis = "x", we fill the array with the corresponding values.
  elseif operator == "Sx"
    for i in 1:n
      m = i - j - 1
      result[i, i+1] = 0.5 * sqrt(j*(j+1) - m*(m+1))
      result[i+1, i] = 0.5 * sqrt(j*(j+1) - m*(m+1))
    end
  # If axis = "y", we fill the array with the corresponding values.
  elseif operator == "iSy"
    for i in 1:n
      m = i - j - 1
      result[i, i+1] = 0.5 * sqrt(j*(j+1) - m*(m+1))
      result[i+1, i] = -0.5 * sqrt(j*(j+1) - m*(m+1))
    end
  # If axis = "z", we fill the array with the corresponding values.
  elseif operator == "Sz"
    for i in 1:n+1
      m = i - j - 1
      result[i, i] = -m
    end
  # Otherwise, we throw an error.
  else
    error("Invalid spin operator")
  end
  return result
end

# With the following code we define the site types and the corresponding spin operators.
# This allows to use the spin operators with higher spin than 1/2 and 1 in the ITensor library.

# We define the site types corresponding to the spin 3/2
ITensors.space(::SiteType"S=3/2") = 4
ITensors.op(::OpName"Sx",::SiteType"S=3/2") = spinOperator("3/2", "Sx")
ITensors.op(::OpName"iSy",::SiteType"S=3/2") = spinOperator("3/2", "iSy")
ITensors.op(::OpName"Sz",::SiteType"S=3/2") = spinOperator("3/2", "Sz")
ITensors.op(::OpName"S+",::SiteType"S=3/2") = spinOperator("3/2", "S+")
ITensors.op(::OpName"S-",::SiteType"S=3/2") = spinOperator("3/2", "S-")

# We define the site types corresponding to the spin 2
ITensors.space(::SiteType"S=2") = 5
ITensors.op(::OpName"Sx",::SiteType"S=2") = spinOperator("2", "Sx")
ITensors.op(::OpName"iSy",::SiteType"S=2") = spinOperator("2", "iSy")
ITensors.op(::OpName"Sz",::SiteType"S=2") = spinOperator("2", "Sz")
ITensors.op(::OpName"S+",::SiteType"S=2") = spinOperator("2", "S+")
ITensors.op(::OpName"S-",::SiteType"S=2") = spinOperator("2", "S-")

# We define the site types corresponding to the spin 5/2
ITensors.space(::SiteType"S=5/2") = 6
ITensors.op(::OpName"Sx",::SiteType"S=5/2") = spinOperator("5/2", "Sx")
ITensors.op(::OpName"iSy",::SiteType"S=5/2") = spinOperator("5/2", "iSy")
ITensors.op(::OpName"Sz",::SiteType"S=5/2") = spinOperator("5/2", "Sz")
ITensors.op(::OpName"S+",::SiteType"S=5/2") = spinOperator("5/2", "S+")
ITensors.op(::OpName"S-",::SiteType"S=5/2") = spinOperator("5/2", "S-")

# We define the site types corresponding to the spin 3
ITensors.space(::SiteType"S=3") = 7
ITensors.op(::OpName"Sx",::SiteType"S=3") = spinOperator("3", "Sx")
ITensors.op(::OpName"iSy",::SiteType"S=3") = spinOperator("3", "iSy")
ITensors.op(::OpName"Sz",::SiteType"S=3") = spinOperator("3", "Sz")
ITensors.op(::OpName"S+",::SiteType"S=3") = spinOperator("3", "S+")
ITensors.op(::OpName"S-",::SiteType"S=3") = spinOperator("3", "S-")

# We define the site types corresponding to the spin 7/2
ITensors.space(::SiteType"S=7/2") = 8
ITensors.op(::OpName"Sx",::SiteType"S=7/2") = spinOperator("7/2", "Sx")
ITensors.op(::OpName"iSy",::SiteType"S=7/2") = spinOperator("7/2", "iSy")
ITensors.op(::OpName"Sz",::SiteType"S=7/2") = spinOperator("7/2", "Sz")
ITensors.op(::OpName"S+",::SiteType"S=7/2") = spinOperator("7/2", "S+")
ITensors.op(::OpName"S-",::SiteType"S=7/2") = spinOperator("7/2", "S-")

# We define the site types corresponding to the spin 4
ITensors.space(::SiteType"S=4") = 9
ITensors.op(::OpName"Sx",::SiteType"S=4") = spinOperator("4", "Sx")
ITensors.op(::OpName"iSy",::SiteType"S=4") = spinOperator("4", "iSy")
ITensors.op(::OpName"Sz",::SiteType"S=4") = spinOperator("4", "Sz")
ITensors.op(::OpName"S+",::SiteType"S=4") = spinOperator("4", "S+")
ITensors.op(::OpName"S-",::SiteType"S=4") = spinOperator("4", "S-")

# We define the site types corresponding to the spin 9/2
ITensors.space(::SiteType"S=9/2") = 10
ITensors.op(::OpName"Sx",::SiteType"S=9/2") = spinOperator("9/2", "Sx")
ITensors.op(::OpName"iSy",::SiteType"S=9/2") = spinOperator("9/2", "iSy")
ITensors.op(::OpName"Sz",::SiteType"S=9/2") = spinOperator("9/2", "Sz")
ITensors.op(::OpName"S+",::SiteType"S=9/2") = spinOperator("9/2", "S+")
ITensors.op(::OpName"S-",::SiteType"S=9/2") = spinOperator("9/2", "S-")

# We define the site types corresponding to the spin 5
ITensors.space(::SiteType"S=5") = 11
ITensors.op(::OpName"Sx",::SiteType"S=5") = spinOperator("5", "Sx")
ITensors.op(::OpName"iSy",::SiteType"S=5") = spinOperator("5", "iSy")
ITensors.op(::OpName"Sz",::SiteType"S=5") = spinOperator("5", "Sz")
ITensors.op(::OpName"S+",::SiteType"S=5") = spinOperator("5", "S+")
ITensors.op(::OpName"S-",::SiteType"S=5") = spinOperator("5", "S-")