# Simple Span Beam Design
By: Margaret Wang

## Design Basis
The program is developed to design a simple span W beam. Gravity loading and axial compression are taken into consideration along the major axis. Compression buckling is also only considered in the major axis. Positive gravity loading is downward, and positive axial load is in compression. Tension loads are not considered, nor are loads applied along the minor axis.

Material properties and design equations are based on W shapes.

## Using the Program
The program is written with Python version 2.7.11, and pandas version 0.17.1

Database of sizes come from AISC: http://www.aisc.org/content.aspx?id=2868

A path to the 'AISC Shapes Database v14.1.csv' should be updated before running the program.

Users are requested for raw input for the area loads, point loads, length, spacing, and effective length factors. Users are also requested to input deflection criteria (default is L/360 and L/240). Inputs should be in integer or float form for loads, length, spacing, and factors. Verification questions should be in character form.

## Methodology
The design process here is based on limiting a database of options based on design requirements, and finding the lightest member that fits that criteria.

The loads are first developed without regard to the added dead load of a beam, which varies by the size chosen. Load and deflection requirements are printed for a basis of comparison. Then, beam weight is added to the loading, and calculated design parameters are added to the database. The final properties can then be associated with each beam shape. 

The method of determining the sizes is a process of elimination. Only shapes that are compliant with flange and web compression are kept. If axial load is present, then the shape ratios for axial compression are included. Only sizes that have a moment of inertia greater than the required moment of inertia for total and live load deflection are kept. The same process repeats for flexural loads, shear loads, and axial loads. The top 10 sizes based on efficiency of section for that loading are listed.

The final result is the lightest shape with the required properties.

## Next Steps
Visualizing the information would be the next step in the process; providing a graphical understanding of where the most likely beams are in terms of strength and serviceability requirements can help an engineer to choose the best solution for the loading conditions.

Additionally, support for minor axis loading, tension requirements, torsional loading, and nonslender/noncompact section results could also be added.

To be a more robust program, support for other shapes should be added.

The next steps could also include developing a program that can allow a person to enter a selected beam shape, and revealing the design performance of that beam, and validating if it works.
