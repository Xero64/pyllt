# Wing Sizing

## Determining the Reference Area

Stall Speed at Maximum Lift Coefficient.

$$
L = W = \frac{1}{2} \cdot \rho \cdot V_{stall}^2 \cdot C_{L,max} \cdot S_{ref}
$$

Therefore for a specific stall speed a specific wing loading must exist.

$$
\frac{W}{S} = \frac{1}{2} \cdot \rho \cdot V_{stall}^2 \cdot C_{L,max}
$$

So therefore $S_{ref}$ can be determined from the mass and stall speed assuming a $C_{L,max}$.

$$
S_{ref} = \frac{2 \cdot W}{\rho \cdot V_{stall}^2 \cdot C_{L,max}}
$$

## Minimising Drag at Cruise

Standard Drag Expression:

$$
D = \frac{1}{2} \cdot \rho \cdot V^2 \cdot C_D \cdot S_{ref}
$$

Drag Coefficient Expression:

$$
C_D = C_{Do} + \frac{{C_L}^2}{\pi \cdot AR \cdot e}
$$

Maximise Lift to Drag Ratio or Minimise Drag to Lift Ratio:

$$
\frac{D}{L} = \frac{C_D}{C_L} = \frac{C_{Do}}{C_L} + \frac{C_L}{\pi \cdot AR \cdot e}
$$

To get $C_L$ at minimum $\frac{D}{L}$, make the deivative equal to zero.

$$
\frac{d(\frac{D}{L})}{d{C_L}} = -\frac{C_{Do}}{{C_L}^2} + \frac{1}{\pi \cdot AR \cdot e} = 0
$$

Therefore:

$$
C_{L,min} = \sqrt{\pi \cdot AR \cdot e \cdot C_{Do}}
$$

$$
C_{D,min} = C_{Do} + \frac{{C_{L,min}}^2}{\pi \cdot AR \cdot e} = 2 \cdot C_{Do}
$$

## Geometry Area at Minimum Drag

So, to determine the aspect ratio for minimum drag at cruise, substitute $C_{L,min}$ as cruise speed and mass.

$$
L_{min} = W = \frac{1}{2} \cdot \rho \cdot V_{cruise}^2 \cdot C_{L,min} \cdot S_{ref} = q_{cruise} \cdot C_{L,min} \cdot S_{ref}
$$

Therefore:

$$
W = q_{cruise} \cdot \sqrt{\pi \cdot AR \cdot e \cdot C_{Do}} \cdot {S_{ref}}
$$

And finally, with re-arranging:

$$
AR = \frac{1}{\pi \cdot e \cdot C_{Do} \cdot q_{cruise}^2} \cdot \left(\frac{W}{S}\right)^2
$$

But aspect ratio is also:

$$
AR = \frac{b^2}{S_{ref}}
$$

Therefore the span $b$ for minimum drag should be:

$$
b = \frac{1}{q_{cruise}} \cdot \sqrt{\frac{S_{ref}}{\pi \cdot e}} \cdot \frac{W}{S}
$$

The mean aerodynamic chord is therefore:

$$
c_{mac} = \frac{S_{ref}}{b}
$$

## Geometry at Minimum Power

For level flight performance:

$$
P = T \cdot V = D \cdot V = \frac{1}{2} \cdot \rho \cdot V^3 \cdot S \cdot \left(C_{Do} + \frac{C_L^2}{\pi \cdot AR \cdot e}\right)
$$

So for 1g flight $L = W$:

$$
C_L = \frac{2 \cdot W}{\rho \cdot V^2 \cdot S}
$$

Therefore power simplifies to:

$$
P = \frac{1}{2} \cdot \rho \cdot V^3 \cdot S \cdot \left(C_{Do} + \frac{4 \cdot W^2}{\pi \cdot AR \cdot e \cdot \rho^2 \cdot V^4 \cdot S^2}\right)
$$

Multiplying out:

$$
P = \frac{1}{2} \cdot \rho \cdot V^3 \cdot S \cdot C_{Do} +  \frac{2 \cdot W^2}{\pi \cdot AR \cdot e \cdot \rho \cdot V \cdot S}
$$

What speed will result in minimum power for a given aspect ratio.

$$
\frac{dP}{dV} = 0 = \frac{3}{2} \cdot \rho \cdot V^2 \cdot S \cdot C_{Do} - \frac{2 \cdot W^2}{\pi \cdot AR \cdot e \cdot \rho \cdot V^2 \cdot S}
$$

Therefore:

$$
V = \sqrt{\sqrt{\frac{4}{3 \cdot \rho^2 \cdot \pi \cdot AR \cdot e \cdot C_{Do} }} \cdot \frac{W}{S}}
$$

Solve for aspect ratio $AR$ given a minimum power air speed $V$:

$$
\pi \cdot AR \cdot e = \frac{1}{3 \cdot q^2 \cdot C_{Do}} \cdot \left(\frac{W}{S}\right)^2
$$