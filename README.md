# Key recovery attacks against the CKKS homomorphic approximate encryption scheme

This repository contains experimental program code implementing our key recovery attacks against the CKKS scheme. Current implementations work with libraries HEAAN, PALISADE, SEAL, HElib, and RNS-HEAAN.

## Build instructions

Makefile expects the HEAAN library to be installed in the current directory.
Install HEAAN and create a symbolic link. For example, if HEAAN is cloned
from github to `/local/HEAAN`, then use this path as `<path-to-HEAAN>`. HEAAN
depends on NTL, which is assumed to have been installed in `/usr/local`

     ln -s <path-to-HEAAN> ./
	 make attack

Similarly, to build the attack for RNS-HEAAN, install RNS-HEAAN (FullRNS-HEAAN)
and create a symbolic link to the current directory. Makefile expects the library
files (.a or .so) are under `./FullRNS-HEAAN/lib`

     ln -s <path-to-FullRNS-HEAAN> ./
	 make rns_attack


For PALISADE, SEAL, and HElib, the Makefile expects their typical installations
into `/usr/local`. Otherwise, you can modify the variable `PALISADE`, `<X>_INCLUDE`,
and `<X>_LIBS` (`X` being `SEAL` or `HELIB`), in Makefile to point to the root include
directory and library directory. These programs also depend on NTL (installed
under `/usr/local` or `/usr`). Then build the executables as:

     make palisade_attack seal_attack helib_attack

Whenever possible, optimized build (e.g. -O3) with parallelization enabled is
preferred for better running times.

## Special note for HElib

The HElib's `Ctxt` class does not provide a public interface for accessing the
component Double-CRT polynomials ("parts" in HElib's terminology). We included
a patch (patches/0001-local-debugging-changes.patch) to add an accessor function `getPart()`
to `Ctxt` and to add a global static variable `decrypted_ptxt_` to store the decrypted
polynomial (for checking encoding error).


## How to run the experiment programs

In general, all programs expect command line arguments to specify the type of
homomorphic computation `<hc>`, the ring dimension `<logN>`, the initial scaling
factor `<logp>`, the upper bound on the random plaintext numbers `<B>`, and the
maximal polynomial degree to evaluate `<deg>`. The homomorphic computation
argument `<hc>` can be one of the following:

     noop, variance, sigmoid, exp

The orders of the arguments are slightly different due to the differences in
how parameters are set up in these libraries. So here are the details:


For HEAAN, run the attack program as:

     ./attack <iter> <hc> <logp> <B> <deg>
     # logN is hard coded in HEAAN's Params.h
     # <iter> is the number of runs to execute, for example, 1

For PALISADE, run palisade_attack as:

     ./palisade_attack <iter> <hc> <logN> <logp> <B> <deg>

For SEAL, run seal_attack as:

     ./seal_attack <iter> <hc> <logN> <logp> <B> <deg>


For HElib, run helib_attack as:

     ./helib_attack <hc> <logN> <logp> <B> <deg>


For RNS-HEAAN, run rns_attack as:

     ./rns_attack <hc> <logN> <L> <logp> <B> <deg>
     # L is the maximal level of computation, for example, 10


For all programs, the parameter `<deg>` is ignored when `<hc>` is noop or variance.
These programs will check and print out the encoding error, and also print out
if the secret key is successfully recovered at the end of a run.
For HEAAN, since power-of-2 modulus is used, with certain probability an inverse
may not exist for the a part of a ciphertext, but the encoding error indicates
if gaussian elimination can be used when collecting a few more ciphertexts.
