
OPTION, PSDUMPFREQ      = 10;
    OPTION, STATDUMPFREQ    = 1;
     OPTION, BOUNDPDESTROYFQ = 10;
    OPTION, AUTOPHASE       = 4;
     
OPTION, VERSION=10900;

Title, string="AWA Photoinjector -- Drive Beamline";

Value,{OPALVERSION};

STRING path_fieldmaps3d="/lstr/sahara/aard/cphillips/awa-opal-lattices/awa_fieldmaps_3D";
STRING path_fieldmaps2d="/lstr/sahara/aard/cphillips/awa-opal-lattices/awa_fieldmaps_2D";

REAL rf_freq             = 1.3e9;
    REAL n_particles         = 3E4;
      REAL beam_bunch_charge   = 1e-9;
     
REAL MX   = 16;
REAL MY   = 16;
REAL MZ   = 32;
REAL BINS =  1;

REAL dTG = 1.0E-11;

REAL Edes    = 1.4e-9;
 REAL gamma   = (Edes+EMASS)/EMASS;
 REAL beta    = sqrt(1-(1/gamma^2));
REAL P0      = gamma*beta*EMASS;
    
value , {Edes, P0};

DiagYAG1: MONITOR, ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0, 3.1245}, OUTFN = "YAG1.h5";
DiagYAG2: MONITOR, ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0, 6.4325}, OUTFN = "YAG2.h5";
DiagYAG3: MONITOR, ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0, 9.6732}, OUTFN = "YAG3.h5";
DiagYAG4: MONITOR, ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,11.5558}, OUTFN = "YAG4.h5";

REAL Gun_field_maps = 2;
  
REAL gun_inj_amp   =  40.64979314068568;
 REAL gun_inj_phase =  -15.985540997564346;
  
photocathode: SOURCE, ORIGIN = {0,0,0};

if (Gun_field_maps == 2.0){   REAL gun_freq = 1300.013338990260;
   REAL gun_len  = 0.292707450;

   GUN: RFCavity, L = gun_len,         VOLT = gun_inj_amp, 	ORIGIN = {0,0,0}, 	TYPE = "STANDING",        FMAPFN = "${path_fieldmaps2d}/driveGun.T7", 	FREQ = gun_freq, 	LAG = (gun_inj_phase*Pi)/180.0;
 }else {   REAL gun_freq = 1300.0;
   REAL gun_len  = 2.327100e-01;

   GUN: RFCavity, L = gun_len,         VOLT = gun_inj_amp, 	ORIGIN = {0,0,0}, 	TYPE = "STANDING",        FMAPFN = "${path_fieldmaps3d}/gun3d.T7", 	FREQ = gun_freq, 	LAG = (gun_inj_phase*Pi)/180.0;
 }
REAL Ifocu = 294.17687528123815;
 REAL KSF = (Ifocu/550.)*0.183100238;
solF: Solenoid, L = 1.2, , KS = KSF,          FMAPFN = "${path_fieldmaps2d}/GUNSOL_F550A.T7";

REAL Imain  = 182.9537919304211;
 REAL KSM = (Imain/440.)*0.610728687;
solM:  Solenoid, L = 1.2, , KS = KSM,          FMAPFN = "${path_fieldmaps2d}/GUNSOL_M440A.T7";

REAL Ibuck = 0.99126080638318*(Ifocu*2.545263e-02+Imain*6.481346e-05)/2.526366e-02;
 REAL KSB = (Ibuck/550.)*0.138950115;
solB:  Solenoid, L = 1.2, , KS = KSB,          FMAPFN = "${path_fieldmaps2d}/GUNSOL_B550A.T7";

REAL Linac_field_maps = 2;
  
REAL Klystron_2 = 1.;
 REAL Klystron_3 = 1.;
 REAL Klystron_4 = 1.;
 
REAL L1_Amp   = 11.779439765020673*Klystron_2;
REAL L1_Phase = -10.885675477233857;
REAL L2_Amp   = 11.186665186098248*Klystron_2;
REAL L2_Phase = 10.59906089386481;
REAL L3_Amp   = 22.06*Klystron_3;
REAL L3_Phase = 0;
REAL L4_Amp   = 22.06*Klystron_4;
REAL L4_Phase = 0;
REAL L5_Amp   = 22.06*Klystron_3;
REAL L5_Phase = 0;
REAL L6_Amp   = 22.06*Klystron_4;
REAL L6_Phase = 0;

REAL Isol_L1 = 0.;
 REAL Isol_L2 = 0.;
 REAL Isol_L3 = 0.;
 
REAL zc_linac1  =  1.1810;
 REAL zc_linac2  =  3.9745;
 REAL zc_linac3  =  5.5425;
 REAL zc_linac4  =  7.6809;
 REAL zc_linac5  =  9.0207;
 REAL zc_linac6  = 10.4358;
 
if (Linac_field_maps == 2.0){    REAL linac_len  = 1.20713;
      REAL linac_freq = 1300.0;

    L1: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",	        FMAPFN = "${path_fieldmaps2d}/driveLinac2D.T7",	VOLT = L1_Amp, LAG = (L1_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,zc_linac1-linac_len/2.};
            
    L2: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",	        FMAPFN = "${path_fieldmaps2d}/driveLinac2D.T7",	VOLT = L2_Amp, LAG = (L2_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,zc_linac2-linac_len/2.};
            
    L3: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",	        FMAPFN = "${path_fieldmaps2d}/driveLinac2D.T7",	VOLT = L3_Amp, LAG = (L3_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,zc_linac3-linac_len/2.};
            
    L4: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",	        FMAPFN = "${path_fieldmaps2d}/driveLinac2D.T7",	VOLT = L4_Amp, LAG = (L4_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,zc_linac4-linac_len/2.};
            
    L5: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",	        FMAPFN = "${path_fieldmaps2d}/driveLinac2D.T7",	VOLT = L5_Amp, LAG = (L5_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,zc_linac5-linac_len/2.};
            
    L6: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",	        FMAPFN = "${path_fieldmaps2d}/driveLinac2D.T7",	VOLT = L6_Amp, LAG = (L6_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,zc_linac6-linac_len/2.};
            }else {     REAL linac_len  = 1.1;
    REAL linac_freq = 1300.400;

    L1: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",		FMAPFN = "${path_fieldmaps3d}/drivelinaccplr_posY.ascii",  	VOLT = L1_Amp, LAG = (L1_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0,-90.0}, ORIGIN = {0,0,zc_linac1-linac_len/2.};
    
    L2: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",		FMAPFN = "${path_fieldmaps3d}/drivelinaccplr_posY.ascii",  	VOLT = L2_Amp, LAG = (L2_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0, 90.0}, ORIGIN = {0,0,zc_linac2-linac_len/2.};
    
    L3: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",		FMAPFN = "${path_fieldmaps3d}/drivelinaccplr_posY.ascii",  	VOLT = L3_Amp, LAG = (L3_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0, 90.0}, ORIGIN = {0,0,zc_linac3-linac_len/2.};
    
    L4: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",		FMAPFN = "${path_fieldmaps3d}/drivelinaccplr_posY.ascii",  	VOLT = L4_Amp, LAG = (L4_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0,-90.0}, ORIGIN = {0,0,zc_linac4-linac_len/2.};
    
    L5: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",		FMAPFN = "${path_fieldmaps3d}/drivelinaccplr_posY.ascii",  	VOLT = L5_Amp, LAG = (L5_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0,-90.0}, ORIGIN = {0,0,zc_linac5-linac_len/2.};
    
    L6: RFCavity, L = linac_len, FREQ = linac_freq, TYPE = "STANDING",		FMAPFN = "${path_fieldmaps3d}/drivelinaccplr_posY.ascii",  	VOLT = L6_Amp, LAG = (L6_Phase* Pi) / 180.0, 	ORIENTATION = {0.0,0.0,-90.0}, ORIGIN = {0,0,zc_linac6-linac_len/2.};
    }
REAL zc_sol_alin1  =  2.0775;
 REAL zc_sol_alin2  =  4.7105;
 REAL zc_sol_alin3  =  6.7772;
 
REAL sol_linac_len = 1.0;
  
REAL KSL1 = (Isol_L1/500.);
REAL KSL2 = (Isol_L2/500.);
REAL KSL3 = (Isol_L3/500.);

sol_L1:  Solenoid, L = 1.0, Z = 0.577435, KS = KSL1,          FMAPFN = "${path_fieldmaps2d}/LINSOL_LB_500A.T7",         ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,zc_sol_alin1-sol_linac_len/2.};
            
sol_L2:  Solenoid, L = 1.0, Z = 3.370935, KS = KSL2,          FMAPFN = "${path_fieldmaps2d}/LINSOL_500A.T7",         ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,zc_sol_alin2-sol_linac_len/2.};
            
sol_L3:  Solenoid, L = 1.0, Z = 4.938935, KS = KSL3,          FMAPFN = "${path_fieldmaps2d}/LINSOL_500A.T7",         ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,zc_sol_alin3-sol_linac_len/2.};
            
DiagYAG0a: MONITOR, ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,0.01}, OUTFN = "YAG0a.h5";
DiagYAG0b: MONITOR, ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,0.6}, OUTFN = "YAG0b.h5";
DiagYAG0c: MONITOR, ORIENTATION = {0.0,0.0,0.0}, ORIGIN = {0,0,zc_sol_alin1-sol_linac_len/2.}, OUTFN = "YAG0c.h5";

DR1: DRIFT, L = 1.4, Z = 0.2;
DR2: DRIFT, L = 3.0, ELEMEDGE = 1.4;

Injector:   Line = (photocathode, GUN, DiagYAG0a, solB, solF, solM, DR1, DiagYAG0b);
Linac:      Line = (L1, DiagYAG0c, sol_L1, DiagYAG1, L2, sol_L2, L3, DiagYAG2,                     sol_L3, L4, L5, DiagYAG3, L6);
 
GSL:        Line = (Injector, L1, DR2);

Gun_to_YAG4:Line = (Injector, Linac, DiagYAG4, 5*DR1);

REAL SigmaLaser = 0.002291188605174065;

DistFT: DISTRIBUTION, TYPE = FLATTOP,        SIGMAX = SigmaLaser,          SIGMAY = SigmaLaser,  
        TRISE = 300e-15/2.35*1.6869,               TFALL = 300e-15/2.35*1.6869,         OFFSETT = 7.000000000000001e-12,        TPULSEFWHM = 6e-12,         CUTOFFLONG = 4.0,        NBIN = BINS,        EMISSIONSTEPS = 200,        EMISSIONMODEL = ASTRA,        EKIN = 0.55,                   EMITTED = True,                WRITETOFILE = True;
    
FS_SC: Fieldsolver, FSTYPE = FFT,             MX = MX, MY = MY, MT = MZ,             PARFFTX = true,             PARFFTY = true,             PARFFTT = true,             BCFFTX = open,             BCFFTY = open,             BCFFTT = open,            BBOXINCR = 1,             GREENSF = INTEGRATED;

BEAM1:  BEAM, PARTICLE = ELECTRON, pc = P0, NPART = n_particles,        BFREQ = rf_freq * 1E-6, BCURRENT = beam_bunch_charge * rf_freq , CHARGE = -1;

TRACK, LINE = Gun_to_YAG4, BEAM = BEAM1, MAXSTEPS = 1000000, DT = {dTg}, ZSTOP={12.0};
 
RUN, METHOD = "PARALLEL-T", BEAM = BEAM1,         FIELDSOLVER = FS_SC, DISTRIBUTION = DistFT;

ENDTRACK;


