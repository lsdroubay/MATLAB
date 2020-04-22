**Lamina Class used in Classical Lamination Theory**

**Author:** Landon Droubay

**Language:** MATLAB


**Description/Purpose:** Returns the needed properties and values for composite material laminate. Uses created [LaminaClass](/LaminaClass) class.

**Output:** This displays the material properties, laminate properties, Qb matrix, loads and moments, 
            ABD matrix, strains and curvatures, normal and shear stresses, normal and shear strain, 
            shear and strain in   , average/smear properties, and the failure factors/criteria.

**Usage/Example:**

For a symmetric carbon fiber 5-layer laminate of angles: 30,-30,0,-30,30. Layer thickness is 1.5e-4 m

```
LaminateAngles = [30,-30,0,-30,30];
symmetric = true;
```
