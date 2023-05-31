

     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 This program allows the csf mixing coefficient and orbital expansion coefficient
 optimization using the graphical unitary group approach and the exponential
 operator mcscf method.
 references:  r. shepard and j. simons, ' int. j. quantum chem. symp. 14, 211 (1980).
              r. shepard, i. shavitt, and j. simons, j. chem. phys. 76, 543 (1982).
              r. shepard in "ab initio methods in quantum chemistry ii" advances in chemical
                  physics 69, edited by k. p. lawley (wiley, new york, 1987) pp. 63-200.
 Original autor: Ron Shepard, ANL
 Later revisions: Michal Dallos, University Vienna

 This Version of Program MCSCF is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.4.0.2     **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 Workspace allocation information:
       222822400 of real*8 words ( 1700.00 MB) of work space has been allocated.

 user input information:

 ======== echo of the mcscf input ========
 ------------------------------------------------------------------------
  &input
   niter=100,
   nmiter=50,
   nciitr=300,
   tol(3)=1.e-5,
   tol(2)=1.e-5,
   tol(1)=1.e-8,
   NSTATE=0,
   npath=1,3,9,10,13,17,19,21,-11,12, 2,
   ncoupl=5,
   tol(9)=1.e-3,
   FCIORB=  1,7,20,1,8,20,1,9,20
   NAVST(1) = 3,
   WAVST(1,1)=1 ,
   WAVST(1,2)=1 ,
   WAVST(1,3)=1 ,
  &end
 ------------------------------------------------------------------------


 ***  Integral file informations  ***


 input integral file : /projects/academic/cyberwksp21/Students/Columbus_tutorial
 /TU

 Integral file header information:
 Hermit Integral Program : SIFS version  srv-p22-12.cbls.c 17:14:25.345 22-Jun-21

 Core type energy values:
 energy( 1)=  3.858870424817E+01, ietype=   -1,    core energy of type: Nuc.Rep.
 total ao core energy =   38.588704248


   ******  Basis set information:  ******

 Number of irreps:                  1
 Total number of basis functions:  36

 irrep no.              1
 irrep label           A  
 no. of bas.fcions.    36


 ***  MCSCF optimization procedure parmeters:  ***


 maximum number of mcscf iterations:        niter=   100

 maximum number of psci micro-iterations:   nmiter=   50
 maximum r,s subspace dimension allowed:    nvrsmx=   30

 tol(1)=  1.0000E-08. . . . delta-emc convergence criterion.
 tol(2)=  1.0000E-05. . . . wnorm convergence criterion.
 tol(3)=  1.0000E-05. . . . knorm convergence criterion.
 tol(4)=  1.0000E-08. . . . apxde convergence criterion.
 tol(5)=  1.0000E-04. . . . small diagonal matrix element tolerance.
 tol(6)=  1.0000E-06. . . . minimum ci-psci residual norm.
 tol(7)=  1.0000E-05. . . . maximum ci-psci residual norm.
 tol(8)=  1.0000E+00. . . . maximum abs(k(xy)) allowed.
 tol(9)=  1.0000E-03. . . . wnorm coupling tolerance.
 tol(10)= 0.0000E+00. . . . maximum psci emergency shift parameter.
 tol(11)= 0.0000E+00. . . . minimum psci emergency shift parameter.
 tol(12)= 0.0000E+00. . . . increment of psci emergency shift parameter.


 *** State averaging informations: ***


 MCSCF calculation performed for  1 DRT.

 DRT  first state   no.of aver.states   weights
  1   ground state          3             0.333 0.333 0.333

 The number of hmc(*) eigenvalues and eigenvectors calculated each iteration per DRT:
 DRT.   no.of eigenv.(=ncol)
    1        4

 orbital coefficients are optimized for the ground state (nstate=0).

 Orbitals included in invariant subspaces:
   symmetry   orbital   mask
       1       7(  7)    20
       1       8(  8)    20
       1       9(  9)    20

 npath(*) options:
  2:  orbital-state coupling terms will be included beginning on iteration ncoupl=  5
  3:  print intermediate timing information.
  9:  suppress the drt listing.
 10:  suppress the hmc(*) eigenvector listing.
 12:  diagonalize the hmc(*) matrix iteratively.
        nunitv= 1 nciitr=** mxvadd=20 nvcimx=20
       rtolci(*),wnorm=     1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 0.0000E+00
   noldv =   0
 13:  get initial orbitals from the formatted file, mocoef.
 17:  print the final natural orbitals and occupations.
 19:  transform the virtual orbitals to diagonalize qvv(*).
 21:  write out the one- and two- electron density for further use (files:mcd1fl, mcd2fl).


   ******  DRT info section  ******


 Informations for the DRT no.  1

 DRT file header:
  title                                                                          
 Molecular symmetry group:    a  
 Total number of electrons:   16
 Spin multiplicity:            1
 Number of active orbitals:    3
 Number of active electrons:   4
 Total number of CSFs:         6
 

 faar:   0 active-active rotations allowed out of:   3 possible.


 Number of active-double rotations:        18
 Number of active-active rotations:         0
 Number of double-virtual rotations:      162
 Number of active-virtual rotations:       81
 lenbfsdef=                131071  lenbfs=                   729
  number of integrals per class 1:11 (cf adda 
 class  1 (pq|rs):         #          21
 class  2 (pq|ri):         #         108
 class  3 (pq|ia):         #         972
 class  4 (pi|qa):         #        1458
 class  5 (pq|ra):         #         486
 class  6 (pq|ij)/(pi|qj): #         342
 class  7 (pq|ab):         #        2268
 class  8 (pa|qb):         #        4374
 class  9 p(bp,ai)         #       13122
 class 10p(ai,jp):        #        2916
 class 11p(ai,bj):        #       15309

 Size of orbital-Hessian matrix B:                    37395
 Size of the orbital-state Hessian matrix C:           4698
 Total size of the state Hessian matrix M:                0
 Size of HESSIAN-matrix for quadratic conv.:          42093


 Source of the initial MO coeficients:

 Input MO coefficient file: /projects/academic/cyberwksp21/Students/Columbus_tutorial/TU
 

               starting mcscf iteration...   1

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58618 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

 trial vectors are generated internally.

 trial vector  1 is unit matrix column     1
 ciiter=   4 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3985323836     -132.9872366318        0.0000000000        0.0000010000
    2       -94.0535597857     -132.6422640338        0.0000000000        0.0000010000
    3       -94.0138933507     -132.6025975989        0.0000000000        0.0000010000
    4       -93.7714030976     -132.3601073458        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  4.288995068762696E-008
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 3.1808E-07 rpnorm= 0.0000E+00 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.790417  -23.169776   -3.089897   -2.348550   -2.033817   -1.932260

 qvv(*) eigenvalues. symmetry block  1
    -0.035872    0.040176    0.084468    0.266514    0.404473    0.997535    1.011503    1.032415    1.185269    1.409150
     1.414620    1.597498    1.808405    1.834155    1.856607    2.487894    2.689971    2.836524    2.931261    3.610651
     3.787297    4.111529    4.414205    4.771371    4.875645    5.505858    5.585668

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    1 emc=    -94.1553285066 demc= 9.4155E+01 wnorm= 3.4312E-07 knorm= 1.0790E-08 apxde= 2.0056E-15    *not conv.*     

               starting mcscf iteration...   2

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 222816479

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 222647099
 address segment size,           sizesg = 222495613
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58618 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   222678332 available sort2 space, avcisx=   222678584

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3985323836     -132.9872366318        0.0000000005        0.0000010000
    2       -94.0535597859     -132.6422640341        0.0000000000        0.0000010000
    3       -94.0138933504     -132.6025975986        0.0000000011        0.0000010000
    4       -93.7714030970     -132.3601073452        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.975772481138430E-008
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 2.4736E-07 rpnorm= 0.0000E+00 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.790417  -23.169776   -3.089897   -2.348550   -2.033817   -1.932260

 qvv(*) eigenvalues. symmetry block  1
    -0.035872    0.040176    0.084468    0.266514    0.404473    0.997535    1.011503    1.032415    1.185269    1.409150
     1.414620    1.597498    1.808405    1.834155    1.856607    2.487894    2.689971    2.836524    2.931261    3.610651
     3.787297    4.111529    4.414205    4.771371    4.875645    5.505858    5.585668

 restrt: restart information saved on the restart file (unit= 13).

 all mcscf convergence criteria are satisfied.

 final mcscf convergence values:
 iter=    2 emc=    -94.1553285066 demc= 5.6843E-14 wnorm= 3.1806E-07 knorm= 8.0817E-09 apxde=-1.1819E-15    *converged*     




   ---------Individual total energies for all states:----------
   DRT #1 state # 1 wt 0.333 total energy=      -94.398532384, rel. (eV)=   0.000000
   DRT #1 state # 2 wt 0.333 total energy=      -94.053559786, rel. (eV)=   9.387186
   DRT #1 state # 3 wt 0.333 total energy=      -94.013893350, rel. (eV)=  10.466565
   ------------------------------------------------------------


 MO-coefficient print-out skipped (no flag 32)
 They may be found in the MOCOEF directory.

          natural orbitals of the final iteration,block  1    -  A  
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.66575125     1.63269253
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.70155623     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30        MO   31        MO   32
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   33        MO   34        MO   35        MO   36
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000
 d1(*), fmc(*), and qmc(*) written to the 1-particle density matrix file.
        105 d2(*) elements written to the 2-particle density matrix file: mcd2fl                                                      


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
    N1_ s       1.998472  -0.000030   1.398408   0.213111   0.000000  -0.000604
    N1_ p       0.000000  -0.000292   0.041269   0.494174   1.315011   0.644431
    N1_ d       0.000000  -0.000105   0.004674   0.005740   0.003805   0.004039
    C1_ s       0.000539   1.999452   0.324541   0.807084   0.000000   0.160542
    C1_ p       0.000476   0.000003   0.098838   0.034165   0.144495   0.703404
    C1_ d      -0.000081   0.000001   0.010188   0.003498   0.003358   0.014052
    H1_ s       0.000301   0.000003   0.052064   0.096151   0.240417   0.066513
    H2_ s       0.000301   0.000003   0.052065   0.096146   0.240418   0.066518
    H3_ s      -0.000004   0.000483   0.008976   0.124965   0.026250   0.170553
    H4_ s      -0.000004   0.000483   0.008977   0.124964   0.026246   0.170552
 
   ao class       7A         8A         9A        10A        11A        12A  
    N1_ s       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
    N1_ p       0.074543   1.424815   0.073979   0.000000   0.000000   0.000000
    N1_ d       0.011182   0.002547   0.003589   0.000000   0.000000   0.000000
    C1_ s       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
    C1_ p       0.904777   0.192226   0.621845   0.000000   0.000000   0.000000
    C1_ d       0.014492   0.013105   0.002143   0.000000   0.000000   0.000000
    H1_ s       0.041375   0.000000   0.000000   0.000000   0.000000   0.000000
    H2_ s       0.041375   0.000000   0.000000   0.000000   0.000000   0.000000
    H3_ s       0.289000   0.000000   0.000000   0.000000   0.000000   0.000000
    H4_ s       0.289006   0.000000   0.000000   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao           N1_        C1_        H1_        H2_        H3_        H4_
      s         3.609356   3.292158   0.496824   0.496826   0.620224   0.620225
      p         4.067930   2.700229   0.000000   0.000000   0.000000   0.000000
      d         0.035473   0.060756   0.000000   0.000000   0.000000   0.000000
    total       7.712759   6.053142   0.496824   0.496826   0.620224   0.620225
 

 Total number of electrons:   16.00000000

 !timer: mcscf                           cpu_time=     0.121 walltime=     0.356
 *** cpu_time / walltime =      0.339
 bummer (warning):timer: cpu_time << walltime.  If possible, increase core memory. event =1
