"""
CS152 - Computational Thinking: Science
Nico Flota Sanchez
08/December/2023

Project_10

This file contains class objects Spiral Galaxy, Star, Simulate, and GUI class objects EntryBox and Button.
"""

import numpy as np
from graphicsPlus import *
import random


class Spiral_Galaxy: 
    def __init__(self, galmass, ahalo, vhalo, rthalo, galpos, galvel): 
        """Defines initial parameters of a spiral galaxy."""
        
        #Dark matter halo
        self.ahalo = ahalo #acceleration
        self.vhalo = vhalo #velocity
        self.rthalo = rthalo #radius 
        
        #Galaxy
        self.galmass = galmass #galaxy mass
        self.galpos = galpos #3D position 
        self.galvel = galvel #velocity
        self.galacc = np.full ((3,1) ,0.) #acceleration 
        #form for array is "(rows, columns), contents".

    def setPosVel(self, position, vel):
        """Sets initial position and velocity for the galaxy"""
        self.galpos = position
        self.galvel = vel
    
    def scaleMass(self, massFact):
        """Scales the overall size of a galaxy as the ratio of the mass of the companion galaxy to the main galaxy. If parameter is equal to 1, then the two galaxies are identical."""
        self.galmass = self.galmass * massFact
        self.vhalo = 1.0 * massFact**0.25  
        self.ahalo = 0.1 * massFact ** 0.5 
        
        a2 = -self.galmass/(self.vhalo**2) 
        a1 = -2 * self.ahalo * self.galmass/(self.vhalo**2) 
        a0 = - self.galmass * (self.ahalo**2)/(self.vhalo**2)
        q = a1 / 3.0 -  (a2**2)/9.0 
        r = (a1 * a2 -3 * a0)/6.0 - (a2**3)/27.0 
        
        gal_1= (r + np.sqrt((q**3)+(r**2)))**0.333 
        gal_2= ((r - np.sqrt((q**3)+(r**2)))**0.333) 
        
        self.rthalo = (gal_1 + gal_2) - a2/3
    
    def MoveGalaxy(self, dt):
        """Evolves the position and velocity during time step 'dt'."""

        new_position = self.galpos + (self.galvel * dt) + (0.5 * self.galacc * (dt**2)) #same equation as in other zelle exercises
        new_vel = self.galvel + (self.galacc * dt)

        self.galpos = new_position
        self.galvel = new_vel

    def Acceleration(self, posin):
        """Acceleration function takes the position of the object (star or galaxy) as input and returns the gravitational acceleration of the object (see Newton's 'g' formula)."""
        #posin = "position of the interacting object" (companion galaxy)

        G = 1.0
        dpos = posin - self.galpos #distance between the two galaxies given as a list of [px, py, pz] values

        r = np.sqrt(np.sum(dpos**2, axis = 0)) #r = np.sqrt((px**2) + (py**2) + (pz**2))

        g = -(G * self.InteriorMass(r)) / (r**2) #newton's formula for gravitational acceleration (g)
        calcacc = (dpos * g) / r

        return calcacc
    
    def Potential(self, posin):
        """Calculates the gravitational potential of the object (star or galaxy) based on the input position."""
        #posin = "position of the interacting object" (companion galaxy)

        G = 1.0
        dpos = posin - self.galpos #distance between the two galaxies given as a list of [px, py, pz] values

        r = np.sqrt(np.sum(dpos**2, axis = 0)) #r = np.sqrt((px**2) + (py**2) + (pz**2))
        
        grav_potential = (-G * self.InteriorMass(r)) / r #Newton's grav potential energy formula (U)

        return grav_potential

    def InteriorMass(self, r):
        """Returns the mass of the object (star or galaxy) based on its position in relation to the center of the galaxy."""

        indices = r < self.rthalo

        intmass = np.full(r.shape, 0.)

        if intmass[indices].shape != (0,):
            #If r < rthalo (r < radius of the dark matter halo) then the companion galaxy is confined within the dark matter halo radius of the galaxy, and so the interior mass is calculated with:

            intmass[indices] = ((self.vhalo**2) * (r[indices]**3)) / ((self.ahalo + r[indices])**2)
        
        if intmass[~indices].shape != (0,):
            #Otherwise the interior mass is simply the mass of the galaxy
            intmass[~indices] = self.galmass
        
        return intmass
    
    def Density(self, r):
        """Determines the density based on the input r (position of the companion object in relation to the center of the main galaxy)"""

        r_inner = r * 0.99
        r_outer = r * 1.01

        m_inner = self.InteriorMass(r_inner)
        m_outer = self.InteriorMass(r_outer)

        diff_mass = (m_outer - m_inner) #difference of masses
        volume = (4/3) * np.pi * ((r_outer**3) - (r_inner ** 3))

        density = diff_mass/volume

        return density
    
    def DynFric(self, pmass, ppos, pvel):
        """Used to consider the role of dynamical friction in our model. It calculates the friction based on the mass, position, and velocity of the object. The resultant friction contributes to the overal acceleration of that object."""

        G = 1.0
        lnGamma = 3.0

        dv = pvel - self.galvel #difference in velocity
        v = np.linalg.norm(dr)

        dr = ppos - self.galpos
        r = np.linalg.norm(dr)

        galrho = self.Density(r)
        fricmag = (4.0 * np.pi * G * lnGamma * pmass * galrho * v) / ((1 + v)**3)
        friction = (-dv / v) * fricmag

        return friction

  
class Star(Spiral_Galaxy):
    """Has the Galaxy class as a parent class."""

    def __init__(self, galmass, ahalo, vhalo, rthalo, galpos, galvel, diskSize, galtheta, galphi, n):
        super().__init__(galmass, ahalo, vhalo, rthalo, galpos, galvel) #inherits all the attributes from the galaxy parent class

        self.diskSize = diskSize #radius in which stars are present within the galaxy
        self.galtheta = galtheta #gives orientation of the galaxy: inclination
        self.galphi = galphi #gives orientation of the galaxy: slew in the orbital plane
        self.n = n #number of stars being simulated. These are evenly split between the two galaxies
        
        #Star templates: creates n stars, each with the position and velocity set to zero
        self.starpos = np.full((3, self.n), 0.)
        self.starvel = np.full((3, self.n), 0.)
        self.staracc = np.full((3, self.n), 0.)

    def MoveStars(self, dt):
        """Evolves the position and velocity of a star during time step 'dt'."""

        new_position = self.starpos + (self.starvel * dt) + (0.5 * self.staracc * (dt**2)) #same equation as in other zelle exercises
        new_vel = self.starvel + (self.staracc * dt)

        self.starpos = new_position
        self.starvel = new_vel

        return self.starpos, self.starvel

    def InitStars(self):
        """Initializes the stars in the galaxy. It generates n number of stars and gives them random positions within the disksize of the galaxy. 
        The velocities of the stars are also calculated based on the positions as well as the input angles phi and theta."""

        cosphi = np.cos(self.galphi)
        sinphi = np.sin(self.galphi)
        
        costheta = np.cos(self.galtheta)
        sintheta = np.sin(self.galtheta)

        #populate the galaxy disk with stars in random positions
        for star in range(self.n):
            bad = True

            while bad:
                xtry = self.diskSize * (1. - 2. * np.random.random())
                ytry = self.diskSize * (1. - 2. * np.random.random())
                rtry = np.sqrt(xtry**2 + ytry**2)
                
                if (rtry < self.diskSize):
                    bad = False
            
            ztry = 0.0

            #convert position to spherical coords
            xrot = xtry*cosphi + ytry*sinphi*costheta + ztry*sinphi*sintheta
            yrot = -xtry*sinphi + ytry*cosphi*costheta + ztry*cosphi*sintheta
            zrot = -ytry*sintheta + ztry*costheta
            
            rot = np.array([xrot, yrot, zrot])
            self.starpos[:, star] = rot + self.galpos.reshape(-1) #changes star position
            

            #calculate escape velocity, get x,y components, and convert to spherical coords
            vesc = np.sqrt(self.InteriorMass(rtry) / rtry) #escape velocity of stars
            
            vxtry = -vesc*ytry / rtry
            vytry = vesc*xtry / rtry
            vztry = 0.0

            vxrot = vxtry*cosphi + vytry*sinphi*costheta + vztry*sinphi*sintheta
            vyrot = -vxtry*sinphi + vytry*cosphi*costheta + vztry*cosphi*sintheta
            vzrot = -vytry*sintheta + vztry*costheta

            vrot = np.array([vxrot, vyrot, vzrot])
            self.starvel[:, star] = vrot + self.galvel.reshape(-1)

            #Acceleration is simply set to zero for the time being
            self.staracc = np.full((1,3), 0.)

        return self.starpos, self.starvel, self.staracc

    def scaleMass(self, massFact):
        """Calculates the disk size of the companion galaxy based of the mass ratio of the companion to the main galaxy."""
        self.diskSize = self.diskSize * np.sqrt(massFact)

        return super().scaleMass(massFact)
    

class Orbit:
    """Orbit class calculates the initial position and velocities of the two galaxies based on the input masses."""

    def __init__(self, energy, rp, tp, eccentricity, mass1, mass2, bod1_pos, bod2_pos, bod1_vel, bod2_vel):
        self.energy = energy
        self.rp = rp #radial distance of closest approach
        self.tp = tp #tangential distance of closest approach
        self.eccentricity = eccentricity
        self.mass1 = mass1
        self.mass2 = mass2
        self.bod1_pos = bod1_pos
        self.bod2_pos = bod2_pos
        self.bod1_vel = bod1_vel
        self.bod2_vel = bod2_vel

        self.initOrbit()
    
    def initOrbit(self):
        #Initializes the orbit of the two galaxies based on the input masses. 
        #This function assumes a parabolic orbit and returns the position and velocity of both galaxies.

        #Parabolic orbit
        mtot = self.mass1 + self.mass2

        p = 2 * self.rp
        nhat = np.sqrt(mtot / (p**3))
        cots = 3.0 * nhat * self.tp
        s = np.arctan(1.0/cots)
        cottheta = (1./(np.tan(s/2)))**0.3333
        theta = np.arctan(1./cottheta)
        tanfon2 = 2. / np.tan(2. * theta)
        r = (p/2.) * (1 + tanfon2**2)

        vel = np.sqrt(2. * mtot / r)
        sin_sq_phi = p / (2. * r)
        phi = np.arcsin(np.sqrt(sin_sq_phi))
        f = 2. * np.arctan(tanfon2)
        
        xc = -r * np.cos(f)
        yc = r * np.sin(f)
        vxc = vel * np.cos(f + phi)
        vyc = -vel*np.sin(f + phi)
        xcom = self.mass2 * xc / (self.mass1 + self.mass2)
        ycom = self.mass2 * yc/ (self.mass1 + self.mass2)
        vxcom = self.mass2 * vxc / (self.mass1 + self.mass2)
        vycom = self.mass2 * vyc / (self.mass1 + self.mass2)

        self.bod1_pos = np.array([[-xcom + self.bod1_pos[0,0]], [-ycom + self.bod1_pos[1,0]], [0.0]])
        self.bod1_vel = np.array([[-vxcom + self.bod1_vel[0,0]], [-vycom + self.bod1_vel[1,0]], [0.0]])
        self.bod2_pos = np.array([[xc - xcom + self.bod2_pos[0,0]], [yc - ycom + self.bod2_pos[1,0]], [0.0]])
        self.bod2_vel = np.array([[vxc - vxcom + self.bod2_vel[0,0]], [vyc - vycom + self.bod2_vel[1,0]], [0.0]])


class EntryBox:
    def __init__(self, win, x, y, width):
        
        point = Point(x, y)

        self.entry = Entry(point, width)
        #entry is (point where it's centered, width)
        #self.entry.setText("0") #starts the box out at zero

        self.entry.draw(win)

        while True:
            key = win.getKey()

            if key == "Return":
                # Get the text from the entry box when Enter is pressed
                break

    def get_value(self):
        """Recovers what was entered into the entry box."""
        return self.entry.getText()


class Button:
    def __init__(self, win, x, y, width, text):
        
        text = Text(Point(x, y), text)
        text.setTextColor("white")
        
        text.draw(win)

        self.button = Rectangle(Point((x - width/2), (y - 12)), Point((x + width/2), (y + 12)))
        self.button.setOutline("white")

        self.ll = self.button.getP1()  # assume p1 is ll (lower left)
        self.ur = self.button.getP2()  # assume p2 is ur (upper right)

        self.button.draw(win)

    
    def button_action(self, win):
        click = False

        while click == False:
            clickPoint = win.getMouse()

            if self.ll.getX() < clickPoint.getX() < self.ur.getX() and self.ll.getY() < clickPoint.getY() < self.ur.getY():
                click = True
                return click

class Simulate:
    def __init__(self, win, comptheta, compphi, galtheta, galphi, peri, massrat, starnum):
        """Creates a simulation by taking parameters given by the user. It uses MakeGalaxies and MakeOrbit to start off the simulation."""
        #User-given parameters
        self.win = win
        self.comptheta = comptheta
        self.compphi = compphi
        self.galtheta = galtheta
        self.galphi = galphi
        self.peri = peri
        self.massrat = float(massrat)
        self.galn = int(0.5 * int(starnum))
        self.compn = int(0.5 * int(starnum))
        
        #will hold Zelle graphics:
        self.vis_bluestars = [] 
        self.vis_yelstars = [] 
        self.vis_galcenter = []
        self.vis_compcenter = []

        self.MakeGalaxies() #saves the position of galactic centers
        self.MakeOrbit() #creates stars in galaxies, draws them, and draws galactic centers

    def MakeGalaxies(self):
        """Creates two galaxies based on the parameters provided in the UI. Draws galactic centers."""
        #Standard parameters: 
        galmass = 4.8
        ahalo = 0.1
        vhalo = 1.0
        rthalo = 5.0
        galpos = np.full((3,1),0)
        galvel = np.full((3,1),25)
        comppos = np.full((3,1),0)
        compvel = np.full((3,1),25)
        diskSize = 50

        #Update position
        galpos[0,0] = 1100
        galpos[1,0] = 200

        comppos[0,0] = 750
        comppos[1,0] = 500

        print(galpos[0,0], galpos[1,0])

        #create the stars in each galaxy
        self.gal = Star(galmass, ahalo, vhalo, rthalo, galpos, galvel, diskSize, self.galtheta, self.galphi, self.galn)
        self.comp = Star(galmass, ahalo, vhalo, rthalo, comppos, compvel, diskSize, self.comptheta, self.compphi, self.compn)

        self.gal.scaleMass(self.massrat) #scales up the main galaxy
      
    def MakeOrbit(self):
        """Creates the parabolic orbit for both galaxies. It then used this orbit to determine the initial position and velocity of the galaxy centers. Finally, it initializes the stars inside both galaxies and then plots the galaxies along with their orbiting stars."""

        energy = 0
        eccentricity = 1
        rperi = 3.0 #radial distance of closest approach

        #given by the user
        tperi = float(self.peri) #tangential distance of closest approach

        print("galpos:", self.gal.galpos, "comppos:", self.comp.galpos)
              
        #Returns the initial position and velocities of our galaxy centers
        self.crashOrbit = Orbit(energy, rperi, tperi, eccentricity, self.gal.galmass, self.comp.galmass, self.gal.galpos, self.comp.galpos, self.gal.galvel, self.comp.galvel)

        print("bod1_pos", self.crashOrbit.bod1_pos, "bod1_vel =", self.crashOrbit.bod1_vel)
        print("bod2,pos", self.crashOrbit.bod2_pos, "bod1_vel =", self.crashOrbit.bod2_vel)

        #Sets the initial position and velocity of the galaxy centers
        self.gal.setPosVel(self.crashOrbit.bod1_pos, self.crashOrbit.bod1_vel) #this sets self.gal.galpos = self.crashOrbit.bod1_pos 
        self.comp.setPosVel(self.crashOrbit.bod2_pos, self.crashOrbit.bod2_vel)

        #Centering the blue galaxy so that the blue (main) galaxy is always in the center of our plot
        """self.graphicsView.setXRange(self.gal.galpos[0,0] - 20, self.gal.galpos[0,0] + 20)
        self.graphicsView.setYRange(self.gal.galpos[1,0] - 20, self.gal.galpos[1,0] + 20)"""

        #Initializes the parameters for the stars in both galaxies
        self.bluestarpos, self.bluevels, self.blueaccs = self.gal.InitStars()
        self.yelstarpos, self.yelvels, self.yelaccs = self.comp.InitStars()

        print(self.bluestarpos[:,0])

        for star in range(self.galn):
            x = self.bluestarpos[0, star]
            y = self.bluestarpos[1, star]
            bluestar = Circle(Point(x, y), 1)
            bluestar.setFill('blue')
            bluestar.setOutline('blue')
            bluestar.draw(self.win)
            self.vis_bluestars.append(bluestar)

        for star in range(self.compn):
            x = self.yelstarpos[0, star]
            y = self.yelstarpos[1, star]
            yelstar = Circle(Point(x, y), 1)
            yelstar.setFill('yellow')
            yelstar.setOutline('yellow')
            yelstar.draw(self.win)
            self.vis_yelstars.append(yelstar)

        #Draw the centers
        galcenter = Circle(Point(self.gal.galpos[0,0], self.gal.galpos[1,0]), 3) 
        #this is Circle(Point([first row, first column], [second row, first column]), radius). thus it gets the x and the y for the Point parameters
        galcenter.setFill("white")
        galcenter.setOutline('white')
        galcenter.draw(self.win)
        self.vis_galcenter.append(galcenter) #add it to the Zelle list of objects

        compcenter = Circle(Point(self.comp.galpos[0,0], self.comp.galpos[1,0]), 3)
        compcenter.setFill("white")
        compcenter.setOutline('white')
        compcenter.draw(self.win)
        self.vis_compcenter.append(compcenter) #add it to the Zelle list of objects

        #dist is the galaxy separation 
        #dist = 3.5 * np.linalg.norm((self.gal.galpos - self.comp.galpos))
    
        #vel is the relative velocity of the two galaxies
        #vel = 250 * np.linalg.norm((self.gal.galvel - self.comp.galvel))                

    def update_plot(self, dt):
        """The update plot function runs the simulation. 
        It calculates the acceleration of the galactic centers for both galaxies as well as the stars within the galaxies. 
        These are then used to move the galaxy centers as well as the stars. The plot is then updated to show this evolution."""

        #We set the restartbutton and pausebutton to True here, while these are active, the plot is constantly updated
        #if the pause button is pressed, the update plot function stops running but we can still see our plots
        #if the reset button is pressed, then we clear all our graphs and reset the time back to zero
        restartbutton = True
        pausebutton = True

        #We always add friction (later we can add functionality to whether or not we want it)
        friction = True

        #dtime taken from Chris Mihos' version, used to evolve the system over time

        if restartbutton == True:
            while pausebutton == True:
                #Calculating the acceleration using our acceleration function for both galaxies. Depends on how far they are from each other
                self.gal.galacc = self.comp.Acceleration(self.gal.galpos)
                self.comp.galacc = self.gal.Acceleration(self.comp.galpos)

                #dist is the relative distance between the two galaxies
                #dist = 3.5 * np.linalg.norm((self.gal.galpos - self.comp.galpos))

                #If Dynamic Friction is present (given by user), then acceleration is calculated slightly differently
                """if friction == True:
                    self.gal.galacc = self.gal.galacc + self.comp.DynFric(self.gal.InteriorMass(dist/3.5), self.gal.galpos, self.gal.galvel)
                    self.comp.galacc = self.comp.galacc + self.gal.DynFric(self.comp.InteriorMass(dist/3.5), self.comp.galpos, self.comp.galvel)
                

                #Acceleration of the center of mass alters each galaxies' accelerations
                comacc = ((self.gal.galmass * self.gal.galacc) + (self.comp.galmass * self.comp.galacc)) / (self.gal.galmass + self.comp.galmass)

                self.gal.galacc = self.gal.galacc - comacc
                self.comp.galacc = self.comp.galacc - comacc"""

                #Calculate acceleration for the stars in galaxy
                self.gal.staracc = self.gal.Acceleration(self.gal.starpos) + self.comp.Acceleration(self.gal.starpos)
                self.comp.staracc = self.comp.Acceleration(self.comp.starpos) + self.gal.Acceleration(self.comp.starpos)     

                #print("galacc:", self.gal.galacc)
                #print("galvel:", self.gal.galvel)

                #After the acceleration is calculated, we then evolve the system using the functions MoveGalaxy and MoveStars
                #print("galpos:", self.gal.galpos)
                self.gal.MoveGalaxy(dt)
                self.bluestarpos, self.bluevels = self.gal.MoveStars(dt) 

                self.comp.MoveGalaxy(dt)
                self.yelstarpos, self.yelvels = self.comp.MoveStars(dt)

                #dist and vel are the rel position and velocity between galaxies
                """ dist = 3.5 * np.linalg.norm((self.gal.galpos - self.comp.galpos))
                vel = 250 * np.linalg.norm((self.gal.galvel - self.comp.galvel))"""

                #again, blue centered is checked, then we want the graph to constantly have the blue (main) galaxy at the focus, so we constantly update the X and Y axis limits so that the galaxy appears at the center
                """if self.bluecenterbox.isChecked():
                    self.graphicsView.setXRange(self.gal.galpos[0,0] - 20, self.gal.galpos[0,0] + 20)
                    self.graphicsView.setYRange(self.gal.galpos[1,0] - 20, self.gal.galpos[1,0] + 20)"""
                
                #Need to move the objects in the visualization:
                #print("one bluestarpos:", self.gal.starpos[:, 0])
                #print(self.comp.starpos[:, 0])
                #print(self.gal.galpos[:, 0])
                #print(self.comp.galpos[:, 0])

                for blue in range(len(self.vis_bluestars)):
                    new_x = self.gal.starpos[0, blue]
                    new_y = self.gal.starpos[1, blue]

                    self.vis_bluestars[blue].move(new_x, new_y)
                
                for yel in range(len(self.vis_yelstars)):
                    new_x = self.comp.starpos[0, yel]
                    new_y = self.comp.starpos[1, yel]

                    self.vis_yelstars[yel].move(new_x, new_y)
                
                #Move the galactic centers
                blue_x = self.gal.galpos[0, 0]
                blue_y = self.gal.galpos[1, 0]

                self.vis_galcenter[0].move(blue_x, blue_y)

                yel_x = self.comp.galpos[0, 0]
                yel_y = self.comp.galpos[1, 0]

                self.vis_compcenter[0].move(yel_x, yel_y)
                


                

        
        

        

        


        

