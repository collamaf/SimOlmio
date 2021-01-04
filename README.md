# Simulation of 166Ho experimental setup

## HOW TO RUN:
```
cd build
cmake -DGeant4_DIR=$G4INSTALL ../
make
./exampleB1

```

## Esempi

### Fantoccio Nema, segnale da Olmio dentro la sfera grande
```
./exampleB1 -Sphere 6 -Source 1 -Isotope 1 -NPrim 1000
```

### Fantoccio Nema, segnale da Tecnezio fuori dalla sfera grande
```
./exampleB1 -Sphere 6 -Source 2 -Isotope 2 -NPrim 1000
```

### Fantoccio Bistecca, segnale da Tecnezio fuori dal cilindro grande
```
./exampleB1 -Sphere -6 -Source 2 -Isotope 2 -NPrim 1000
```

### Misure all'IFO per WIDMApp (Dec 2020):
- Sfera r=8.5mm con 99mTc:
```
./exampleB1 -Sphere 3 -Source 1  -SphereDistY 55 -DetConf 1 -Isotope 2 -NPrim 1000000
```

- Sfera r=11mm con 64Cu:
```
./exampleB1 -Sphere 4 -Source 1  -SphereDistY 70 -DetConf 3 -Isotope 4 -NPrim 1000000
```

- Sfera r=1mm con 131I (approssimazione rozza della pasticca di Iodio):
```
./exampleB1 -Sphere 1 -Source 1  -SphereDistY 11 -DetConf 2 -Isotope 3 -NPrim 1000000
```


### Source Choice:
1 - Activate spheres/cylinders
2 - Activate phantom
This is achieved by means of the "GPS/confine" function

### Geometry Choice is performed by means of -Sphere argument:
- -1/6 NEMA-like phantom with increasing diameter sphere in the center
- -1/-6 "Steak"-like phantom with increasing diameter inner cylinder (steak_thickness hard coded to 3cm)



## OUTPUT:
A root file named CMOSmc_{XX}.root is created, reporting the several parameters used for the run, in which on an event (i.e. a primary particle) by event basis it is stored:

### SOURCE vector (one entry per primary particle, only first 100k events are written for disk space sake):
- AllX: X coordinate of primary particle [mm];
- AllY: Y coordinate of primary particle [mm];
- AllZ: Z coordinate of primary particle [mm];
- AllCosX[]: X directive cosine of produced electron;
- AllCosY[]: Y directive cosine of produced electron;
- AllCosZ[]: Z directive cosine of produced electron;
- AllEne[]: kinetic energy of produced particle  [keV];
- AllIsotope[]: parentID - 1 of track: so for example for Sr source is 0 if track is son of Sr, 1 if of Y;
- ExitX[]: X coordinate of primary particle exiting the source volume [mm]; ("Exiting" means going from source to dummy or from absorber to dummy if there is an absorber)
- ExitY[]: Y coordinate of primary particle exiting the source volume [mm];
- ExitZ[]: Z coordinate of primary particle exiting the source volume [mm];
- ExitCosX[]: X directive cosine of primary particle exiting the source volume;
- ExitCosY[]: Y directive cosine of primary particle exiting the source volume;
- ExitCosZ[]: Z directive cosine of primary particle exiting the source volume;
- ExitEne[]: kinetic energy of primary particle exiting the source volume [keV];
- ExitPart[]: kind of primary particle (11=e-, -11=e+, 22=gamma, 13=mu-...) exiting the source volume;
- ExitParentID[]: partent-id of particle exiting the source
- ExitProcess[]: process that created the particles that exits the source (see table above)
- ExitTrackN: number of different tracks exiting the source per event

### B1 vector (one entry per primary particle):
- Eabs: energy absorbed in CMOS [keV];
- EAbsComp[2]: vector containing energy absorbed in CMOS [keV] due to Sr (comp 1) and to Y (comp 2)
- PreGCTrackN: number of tracks entering the GC (if present, empty otherwise);
- PreGCPart: kind of particle entering the GC (if present, empty otherwise);
- PreGCEn: energy of particle entering the GC (if present, empty otherwise) [keV];
- PreGCX: x position of each tracks entering GC [mm];
- PreGCY: y position of each tracks entering GC [mm];
- PreGCZ: z position of each tracks entering GC [mm];
- PreGCDirX: x direction of each tracks entering GC;
- PreGCDirY: y direction of each tracks entering GC;
- PreGCDirZ: z direction of each tracks entering GC;
- PreGCVX: x position of each vertex of tracks entering GC [mm];
- PreGCVY: y position of each vertex of tracks entering GC [mm];
- PreGCVZ: z position of each vertex of tracks entering GC [mm];
- PreGCEnPrim: energy of the primary particle that originated the track that is now entering GC [keV];
- PreGCEventNum: number of times that the event of the primary particle that originated the track that is now entering GC  gave a track entering GC (useful to check for double countings);
- PreGCOrigReg: Number of the region from which the primary whose child now enters the GC originated [4: Phantom Shell, 3: Phantom, 2: Sphere Shell, 1: Sphere];
- PreGCPrestepReg: Number of the preStep region [4: Phantom Shell, 3: Phantom, 2: Sphere Shell, 1: Sphere];
- PreGCCreatorProcess: Creator process of track entering GC [code];
- PreGCCptReg: Region in which happened compton scattering;
- SourceX: X coordinate of primary particle (isotope) giving a signal in GC [mm];
- SourceY: Y coordinate of primary particle (isotope) giving a signal in GC [mm];
- SourceZ: Z coordinate of primary particle (isotope) giving a signal in GC [mm];
- SourceCosX[]: X directive cosine of decay electron(s) giving a signal in GC;
- SourceCosY[]: Y directive cosine of  decay electron(s) giving a signal in GC;
- SourceCosZ[]: Z directive cosine of decay electron(s) giving a signal in GC;
- SourceEne[]: kinetic energy of produced particle  [keV];
- SourcePart: primary particle;
- SourceReg: region of primary particle;


## CHANGELOG
2020.12.23 by collamaf
- UPDATE Readme from old CMOS sim
- Manually merge root files at the end (and delete temporary root files) since the automatic G4 way fails to merge vectors
- New argument: SphereDistY. If !=0 is the distance from the NEMA surface to the sphere center along Y
- Possibility to add detectors like IFO WIDMApp studies: if DetConf!=0 3 detectors are placed, and 3 configs can be chosen via DetConf1, 2 or 3. For now score only particles entering these hollow detectors.
- Add 64Cu and 131I sources

2020.12.29 by collamaf
- Reorganization of geometry, now with 2 different DetConst, one for NEMA and one for Steak (previous selection mechanism kept)
- In Steak Geometry now changed whole orientation (towards Z), renamed volumes and added triangle inner vessel if -Sphere <-6 (for now triangles dimensions are fixed).

2021.01.04 by collamaf
- Add PTER and deposited energy scoring for detectors in WIDMApp-like case

## TO DO's


