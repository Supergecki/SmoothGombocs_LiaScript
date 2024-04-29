# An introduction to smooth Gombocs (by Martin Wille)

## Disclaimer
This presentation is based on a mathematical project I did during 11th and 12th grade. Just in case this presentation goes too much into the mathematical direction, remind me to explain things more thoroughly if you didn't get them right away.

## Motivation
### Tumbler toys
      ![Tumbler toy](https://i.makeagif.com/media/4-26-2023/XcbA8k.gif)

{{1}}

* When leaning the toy to the side, it makes itself stand upright again.
* Witchcraft?!

            {{2}}
* Certainly not, it's Physics!
* Center of mass near the bottom stabilizes the toy into the upright position.
* Realized by material of higher density (e.g. iron) at the bottom.
* What could a mathematician ask next...?

{{3}}
**Can we do this with homogenous density?**

### The Gomboc
Yes! Meet the gomboc:
!?[Gomboc in action](https://www.youtube.com/watch?v=_VuGW_WQG7A)

{{1}}
* Only the surface shape causes it to stabilize itself!
* Witchcraft?!

{{2}}
* Depends on who you ask...
* ...it's interesting nevertheless.

## Problem

What's the problem now?

{{1}}
<section>
![gomboc still](https://gomboc-shop.com/sites/default/files/gomboc.jpg)

See the edges and vertices?
</section>

{{2}}
Mathematicians would rather not have them. 
_(For calculus reasons.)_

{{3}}
=> Introducing _smooth gombocs_ (without edges and vertices).

## The formula

      {{1}}
Luckily, we do have one simple formula for smooth gombocs!

      {{2}}
$$
   R(\theta,\varphi,c,d)=1+d\cdot \Delta R(\theta,\varphi,c)\\
$$

      {{3}}
...did I say one?
$$
   \Delta R(\theta,\varphi,c)=\begin{cases}a\cdot f_1+(1-a)\cdot f_2,& |\varphi|<\pi/2\\
   1,&\varphi=\pi/2\\
   -1,&\varphi=-\pi/2
   \end{cases}
$$

      {{4}}
...I didn't say 'simple', did I?
$$
   a(\theta,\varphi,c)=\dfrac{\cos^2(\theta)\cdot(1-f_1^2)}{\cos^2(\theta)\cdot(1-f_1^2)+\sin^2(\theta)\cdot(1-f_2^2)}
   $$$$
   f_1(\varphi,c)=\sin(f(\varphi,c))$$$$
   f_2(\varphi,c)=-f_1(-\varphi,c)$$$$
   f(\varphi,c)=\pi\cdot\left(\dfrac{e^{\frac{\varphi}{\pi c}+\frac{1}{2c}}-1}{e^{1/c}-1}-\dfrac{1}{2}\right)
$$

### To be fair...

...I don't really know what most of this pile of notation does either.

      {{1}}
Visualization time!

      {{2}}
We'll look at the Python code behind it.

## Let's get into the code

{{1}}
* What will the code do? ==> Calculate and generate data points that can be visualized in 3D space using an external library.

{{2}}
* How does it do that? ==> Most of the code is calculation according to the formulas.

{{3}}
* And the rest? ==> Constants, file handling and debug outputs. Not much more.

{{4}}
<section>
* Can we see it? ==> Of course!

``` python heightfield.py
# Creates a heightfield and fills it with values. Uses the Matrix class.

from matrix import Matrix
from math import pi, sin, cos, tan, exp
import time

RAD = 180           # Represents a half angle, corresponding to pi in radiant.
PARAM_C = 1/2       # smaller values create bigger deviation of base_f from the identity function
PARAM_D = 1/2       # deviation from the sphere
ANGLE_ACCURACY = 180 # The higher the angle accuracy, the more values are calculated. It is the number of phi respectively theta values.
VALUE_ACCURACY = 3  # Accuracy of the values in number of decimals after the decimal point.

PRINT_HEIGHTFIELD = True # If set to True, prints the heightfield matrix into the console.
SAVE_HEIGHTFIELD  = True  # If set to True, saves a copy of the heightfield in a text file.
FILENAME_HF       = "heightfield1.txt" # Filename for saving.
CREATE_XYZ_FILE   = True  # If set to True, creates a xyz-file to plot the point cloud.
FILENAME_XYZ      = "test.xyz" # Filename for saving.
PRINT_TIME        = False # If set to True, prints out used time.

def get_phi(m):
    return -RAD/2 + m * (RAD / ANGLE_ACCURACY)
    
def get_theta(n):
    return n * (2*RAD / ANGLE_ACCURACY)
    
def convert_to_rad(angle):
    return angle * pi / 180
 
def base_f(phi, c):
    return pi * ((exp(phi / (RAD * c) + 1 / (2 * c)) - 1) / (exp(1/c) - 1) - 1/2) 
    
def f1(phi, c):
    return sin(base_f(phi, c))
    
def f2(phi, c):
    return -f1(-phi, c)
    
def calculate_R(phi, theta, c, d):
    return 1 + d * calculate_Rdelta(phi, theta, c)
    
def calculate_Rdelta(phi, theta, c):
    if phi == RAD/2:
        return 1
    elif phi == -RAD/2:
        return -1
    else:
        a = get_weighting(phi, theta, c)
        return a * f1(phi, c) + (1 - a) * f2(phi, c)
        
def get_weighting(phi, theta, c):
    rad_phi = convert_to_rad(phi)
    rad_theta = convert_to_rad(theta)
    return 1 / (1 + tan(rad_theta)**2 * cos(base_f(phi, c))**2 / cos(base_f(-phi, c))**2)
    #return (cos(rad_theta)**2 * (1 - f1(phi, c)**2))   /   (cos(rad_theta)**2 * (1 - f1(phi, c)**2) + sin(rad_theta)**2 * (1 - f2(phi, c)**2))

now = time.time()
heightfield1 = Matrix(ANGLE_ACCURACY, ANGLE_ACCURACY)
if CREATE_XYZ_FILE: point_cloud1 = Matrix(ANGLE_ACCURACY**2, 3)
for m in range(len(heightfield1)):
    for n in range(len(heightfield1[m])):
        radius = calculate_R(get_phi(m), get_theta(n), PARAM_C, PARAM_D)
        heightfield1[m][n] = round(radius, VALUE_ACCURACY)
        if CREATE_XYZ_FILE:
            phi = convert_to_rad(get_phi(m))
            theta = convert_to_rad(get_theta(n))
            #radius = 1
            point_cloud1[m * len(heightfield1) + n][0] = round(radius * cos(theta) * cos(phi), VALUE_ACCURACY) # x coordinate for point
            point_cloud1[m * len(heightfield1) + n][1] = round(radius * sin(theta) * cos(phi), VALUE_ACCURACY) # y coordinate for point
            point_cloud1[m * len(heightfield1) + n][2] = round(radius * sin(phi) , VALUE_ACCURACY)           # z coordinate for point

if PRINT_HEIGHTFIELD:
    print("Values for theta, ranging from 0 to 2*pi:".center(ANGLE_ACCURACY * (VALUE_ACCURACY + 3) + (ANGLE_ACCURACY - 1)))
    print("\n")
    heightfield1.fprint()
    print("\n")
    print("Side values are for phi, ranging from -pi/2 to pi/2.".center(ANGLE_ACCURACY * (VALUE_ACCURACY + 3) + (ANGLE_ACCURACY - 1)))
    
if SAVE_HEIGHTFIELD:
    f = open(FILENAME_HF, 'w')
    f.write(str(heightfield1))
    f.close()
    
if CREATE_XYZ_FILE:
    f = open(FILENAME_XYZ, 'w')
    f.write(str(point_cloud1))
    f.close()

if PRINT_TIME:
    print("Time used (seconds): " + str(time.time()-now))

```

</section>

### Imports and parameters
The program can be modified by changing these parameters.

```python heightfield.py
from matrix import Matrix
from math import pi, sin, cos, tan, exp
import time

RAD = 180           # Represents a half angle, corresponding to pi in radiant.
PARAM_C = 1/2       # smaller values create bigger deviation of base_f from the identity function
PARAM_D = 1/2       # deviation from the sphere
ANGLE_ACCURACY = 180 # The higher the angle accuracy, the more values are calculated. It is the number of phi respectively theta values.
VALUE_ACCURACY = 3  # Accuracy of the values in number of decimals after the decimal point.

PRINT_HEIGHTFIELD = True # If set to True, prints the heightfield matrix into the console.
SAVE_HEIGHTFIELD  = True  # If set to True, saves a copy of the heightfield in a text file.
FILENAME_HF       = "heightfield1.txt" # Filename for saving.
CREATE_XYZ_FILE   = True  # If set to True, creates a xyz-file to plot the point cloud.
FILENAME_XYZ      = "test.xyz" # Filename for saving.
PRINT_TIME        = False # If set to True, prints out used time.
```

### Calculations
The formulas shown earlier and nothing else.

```python heightfield.py
def get_phi(m):
    return -RAD/2 + m * (RAD / ANGLE_ACCURACY)
    
def get_theta(n):
    return n * (2*RAD / ANGLE_ACCURACY)
    
def convert_to_rad(angle):
    return angle * pi / 180
 
def base_f(phi, c):
    return pi * ((exp(phi / (RAD * c) + 1 / (2 * c)) - 1) / (exp(1/c) - 1) - 1/2) 
    
def f1(phi, c):
    return sin(base_f(phi, c))
    
def f2(phi, c):
    return -f1(-phi, c)
    
def calculate_R(phi, theta, c, d):
    return 1 + d * calculate_Rdelta(phi, theta, c)
    
def calculate_Rdelta(phi, theta, c):
    if phi == RAD/2:
        return 1
    elif phi == -RAD/2:
        return -1
    else:
        a = get_weighting(phi, theta, c)
        return a * f1(phi, c) + (1 - a) * f2(phi, c)
        
def get_weighting(phi, theta, c):
    rad_phi = convert_to_rad(phi)
    rad_theta = convert_to_rad(theta)
    return 1 / (1 + tan(rad_theta)**2 * cos(base_f(phi, c))**2 / cos(base_f(-phi, c))**2)
    #return (cos(rad_theta)**2 * (1 - f1(phi, c)**2))   /   (cos(rad_theta)**2 * (1 - f1(phi, c)**2) + sin(rad_theta)**2 * (1 - f2(phi, c)**2))
```

### Calculation and file handling
Prepares the data structures for the calculation and is able to output the results to files or the command line, depending on your parameters.

```python heightfield.py
now = time.time()
heightfield1 = Matrix(ANGLE_ACCURACY, ANGLE_ACCURACY)
if CREATE_XYZ_FILE: point_cloud1 = Matrix(ANGLE_ACCURACY**2, 3)
for m in range(len(heightfield1)):
    for n in range(len(heightfield1[m])):
        radius = calculate_R(get_phi(m), get_theta(n), PARAM_C, PARAM_D)
        heightfield1[m][n] = round(radius, VALUE_ACCURACY)
        if CREATE_XYZ_FILE:
            phi = convert_to_rad(get_phi(m))
            theta = convert_to_rad(get_theta(n))
            point_cloud1[m * len(heightfield1) + n][0] = round(radius * cos(theta) * cos(phi), VALUE_ACCURACY) # x coordinate for point
            point_cloud1[m * len(heightfield1) + n][1] = round(radius * sin(theta) * cos(phi), VALUE_ACCURACY) # y coordinate for point
            point_cloud1[m * len(heightfield1) + n][2] = round(radius * sin(phi) , VALUE_ACCURACY)           # z coordinate for point

if PRINT_HEIGHTFIELD:
    print("Values for theta, ranging from 0 to 2*pi:".center(ANGLE_ACCURACY * (VALUE_ACCURACY + 3) + (ANGLE_ACCURACY - 1)))
    print("\n")
    heightfield1.fprint()
    print("\n")
    print("Side values are for phi, ranging from -pi/2 to pi/2.".center(ANGLE_ACCURACY * (VALUE_ACCURACY + 3) + (ANGLE_ACCURACY - 1)))
    
if SAVE_HEIGHTFIELD:
    f = open(FILENAME_HF, 'w')
    f.write(str(heightfield1))
    f.close()
    
if CREATE_XYZ_FILE:
    f = open(FILENAME_XYZ, 'w')
    f.write(str(point_cloud1))
    f.close()

if PRINT_TIME:
    print("Time used (seconds): " + str(time.time()-now))

```

### ...anything visual?
See program running on my PC... If that does not work, here are some pictures:
![c=1/2 d=1/2 gomboc](https://raw.githubusercontent.com/Supergecki/SmoothGombocs_LiaScript/main/picture%20files/c0%2C5d0%2C5.PNG "Gomboc generated with c=1/2 and d=1/2.")

![c=1/2 d=1 gomboc](https://raw.githubusercontent.com/Supergecki/SmoothGombocs_LiaScript/main/picture%20files/c0%2C5d0%2C5.PNG "Gomboc generated with c=1/2 and d=1.")

![c=2 d=1/2 gomboc](https://raw.githubusercontent.com/Supergecki/SmoothGombocs_LiaScript/main/picture%20files/c0%2C5d1.PNG "Gomboc generated with c=2 and d=1/2.")

## And now?
* Original project includes program to analyze these gomboc structures.
* This is just an abstract of the full "paper". Going too much into detail would take too long!

      {{1}}
**Thank you for your attention and for following me into this little mathematical world!** 

{{2}}
Have fun with this short video of a tortoise behaving like a gomboc...
!?[gomboc tortoise](https://plus.maths.org/issue52/features/gomboc/Gombocturtle.m4v)

## Sources
* Made with LiaScript! Just in case anyone wants to read the documentation: https://liascript.github.io/course/?https://raw.githubusercontent.com/liaScript/docs/master/README.md#1
* Stand-up toy GIF: https://makeagif.com/gif/stehaufmannchen-XcbA8k
* Gomboc video: https://www.youtube.com/watch?v=_VuGW_WQG7A
* Gomboc photo: https://gomboc-shop.com/sites/default/files/gomboc.jpg
* For other sources, see my BeLL "Über die Existenz glatter Gömböcs" (2022) and the sources mentioned there.
* GitHub repository containing the code shown (and some more): https://github.com/Supergecki/smooth_gombocs