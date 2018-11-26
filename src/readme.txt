Modern Fortran Workshop 26.11.13

European project: PoP Coe  help developing parralel codes for free.

Fortran 2003 and then 2008: introduce object oriented prog and parallelism (no lib) CoArrays?
Fortran 2018: no major changes (improve CoArray features)

CUDA fortran, for example, is an extension of the languague (not part of the official language) => No guarantee it will always be supported.

Best practices:

1) label IF and LOOP blocks

LABEL: if() then

endif LABEL

2) capital letters only for Subroutines and Functions and use space between arguments for Sub and Func. (not for arrays) => easier to read

3) Fortran is column-major => A(:,i) faster than A(i,:) => better cache-usage + vectorisation

4) Capitalise names of constants

integer, parameter :: MAX_CELLS = 1000

5) Comments in the beginnig of the file with
purpose
name
date
licence

6) meaninigful names

7) use intent keyword (intent(out) by default

8) capitalize Library functions

9) Fortran IF check all tests eg if( a==1 .and b==2) check both better split into two if sometimes (segmentation fault for example)

10) Use lbound and ubound to get  array bounds (can pass to a subroutine without their bounds (if one uses contiguous keyword) )

11) Use brackets to indicate arrays eg matrix(:,:) better than matrix

12) DO NOT use do loop for array assignments or array operations => avoid bugs and compilers can vectorise

13) Use gdb for debugging (what is gdb?) 

14) Use intrinsic functions like where, all, sim, log, abs, pack, transpose, reshape ... (vectorisation)

15) some compilers like intel have options for vectorisation and optimisation

16) cshift? all?

17) derived types use the t suffix
    type point_t tells you it is a derived type

18) asiign type p1 = point_t(arg1,arg2,...)

19) pointers use p suffix center_p

20) extends a derived type => oop heritage

21) Derived data types can be I/O

22) Fortran 2003 introduced command line arguments  

23) DO NOT USE GOTO => Use cycle or exit in loops

24) Fortran Block statement => organise the code as Blocks

25) DO NOT NEGLECT DOCUMENTATION: every program, module, subroutine, LaTex syntax should be used

26) DO NOT ALLOCATE large arrays in subroutine because they are allocated on the stack. ALLOCATE in the main program to use the "heap". What is HEAP? => heap for global variables

27) Use Pre-Processing

#ifdef DEBUG
  print ....
#endif

nagfor -c -DDEBUG code.F90

28) Preprocessing using Fypp => can use Python modules within the Fortran code

29) use, intrinsic :: iso_fortran_env to use single, double,... always postfix _DP

30) quadruple precision can be employed REAL128  (be careful with GNU Fortran?)

31) Polymorphism => if one has several routines doing the same things but on different types it is possible to give a single same and let the code use the right routine
 
    inside a module:

    interface my_sum
       module procedure real_sum
       module procedure int_sum
    end interface

    now define functions/subroutines real_sum and int_sum

32) Fortran submodules => separate a module into two files => one for the declaration of the subroutines and one for the implementation of the subroutines (no need to recompile + may avoid bugs)

33) Avoid do while loops because they cannot be parallelised. 

34) do concurent loop => all iterations in the loop are independant => allows vectorisations (no always openMP available)

35) optional arguments in Fortran
    real, intent(in), optional :: .... 

36) declare pure and elemental functions and subroutines help compiler optimisations

37) when building libraries, use the name of the libraries as prefix
    eg
    use HAWK
    call HAWK_init ....

38) nagfor can automatically generate modern fortran codes from old ones using polish option

39) link with lib.a => static linking = link only the few subroutines that are used ie not the whole library and faster than shared; with lib.so => shared linking (must define LD_LIBRARY_PATH)

40) ldd executable tells you the shared lib that are required 

41) As a good practice, use 2 compilers one for portability (eg NAG) and one for performance (eg Intel)

42) makedef90 gives you all the dependencies which helps to write a Makefile 

43) CMake is also good mature option 

44) Fortran Documenter (FORD)  => dependencies are visualised (on GitHub), support LaTeX equations, ... => create HTML documentation linked to the source code

45) LaTex Package for Printing Fortran code <=> Listing

