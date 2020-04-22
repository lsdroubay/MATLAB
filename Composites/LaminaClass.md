**Lamina Class used in Classical Lamination Theory**

**Author:** Landon Droubay

**Language:** MATLAB


**Description/Purpose:** Class to hold the needed properties and values for composite material lamina. 

**Usage:** See example of it's use in my [Classical Lamination Theory Analysis](/CLT)

```
classdef LaminaClass
   properties
      material
      theta
      thickness
      E1
      E2
      v12
      v21
      G12
      Q
      S
      a
      axy
      %Failure Stresses
      sigmaT1
      sigmaC1
      sigmaT2
      sigmaC2
      tauF12
   end
   methods
    
      
   end
end
```
