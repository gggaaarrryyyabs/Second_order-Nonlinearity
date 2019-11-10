********
Mini-AES
********

Description
===========

Raphael Chung-Wei Phan presented a version of the AES [Phan02miniadvanced]_, with all the parameters significantly reduced while preserving its original structure. This Mini version is purely educational and is designed to grasp the underlying concepts of Rijndael-like ciphers. It may also serve as a testbed for starting cryptanalysts to experiment with various cryptanalytic attacks. The Mini-AES cipher is a :math:`16 \times 16` vector Boolean function and the Mini-AES encryption is performed with a secret key of 16 bits.

The Mini-AES S-box is called *NibbleSub*, and defines a simple operation that substitutes each input with an output according to a 4x4 substitution table (S-box) given in the table below. These values are, in fact, taken from the first row of the first S-box in DES.

+-----------------------+
| NibbleSub Truth Table |
+=======+===============+
| Input | Output        |
+-------+---------------+
| 0000  | 1110          |
+-------+---------------+
| 0001  | 0100          |
+-------+---------------+
| 0010  | 1101	        |
+-------+---------------+
| 0011  | 0001	        |
+-------+---------------+
| 0100  | 0010	    	|
+-------+---------------+
| 0101  | 1111          |
+-------+---------------+
| 0110  | 1011 	        |
+-------+---------------+
| 0111  | 1000          |
+-------+---------------+ 
| 1000  | 0011	        |
+-------+---------------+
| 1001  | 1010	        |
+-------+---------------+
| 1010  | 0110 	        |
+-------+---------------+
| 1011  | 1100          |
+-------+---------------+
| 1100  | 0101          |
+-------+---------------+
| 1101  | 1001	        |
+-------+---------------+
| 1110  | 0000	        |
+-------+---------------+
| 1111  | 0111	        |
+-------+---------------+

The inverse of the previous table is easily computed by interchanging the input nibble with the output nibble, and then resorting it based on the new input nibble, as given in the table below.

+--------------------------+
| NibbleSubInv Truth Table |
+=======+==================+
| Input | Output           |
+-------+------------------+
| 0000  | 1110             |
+-------+------------------+
| 0001  | 0011             |
+-------+------------------+
| 0010  | 0100             |
+-------+------------------+
| 0011  | 1000             |
+-------+------------------+
| 0100  | 0001             |
+-------+------------------+
| 0101  | 1100             |
+-------+------------------+
| 0110  | 1010             |
+-------+------------------+
| 0111  | 1111             |
+-------+------------------+
| 1000  | 0111             |
+-------+------------------+
| 1001  | 1101             |
+-------+------------------+
| 1010  | 1001             |
+-------+------------------+
| 1011  | 0110             |
+-------+------------------+
| 1100  | 1011             |
+-------+------------------+
| 1101  | 0010             |
+-------+------------------+
| 1110  | 0000             |
+-------+------------------+
| 1111  | 0101             |
+-------+------------------+

Summary
=======

+--------------+------+------+-------+------+-------+------+---------+----------------+--------+------+
| S-box        | size | *NL* | *NL2* | *LD* | *DEG* | *AI* | *MAXAC* | :math:`\sigma` | *LP*   | *DP* |
+==============+======+======+=======+======+=======+======+=========+================+========+======+
| NibbleSub    | 4x4  | 2    | 0     | 0    | 2     | 2    | 16      | 1408           | 0.5625 | 0.5  |
+--------------+------+------+-------+------+-------+------+---------+----------------+--------+------+
| NibbleSubInv | 4x4  | 2    | 0     | 0    | 2     | 2    | 16      | 11408          | 0.5625 | 0.5  |
+--------------+------+------+-------+------+-------+------+---------+----------------+--------+------+
| MixColumn    | 8x8  | 0    | -     | 0    | 1     | 1    | 256     | 16777216       | 1      | 1    |
+--------------+------+------+-------+------+-------+------+---------+----------------+--------+------+

NibbleSub
=========

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSub.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSub.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSub.anf>`_

`Characteristic Function <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSub.char>`_

`Walsh Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSub.wal>`_

+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|16|0 |0 |0 |0 |0 |0 |0 |0 |0  |0 |0 |0 |0 |0  |0  |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |0 |-4|-4|0 |0 |-4|12|4 | 4 |0 |0 |4 |4 |0  |0  |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |0 |-4|-4|0 |0 |-4|-4|0 |0  |4 |4 |0 |0 |-12|4  |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |0 |0 |0 |0 |0 |0 |0 |4 |-12|-4|-4|4 |4 |-4 |-4 |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |4 |0 |-4|-4|-8|-4|0 |0 |-4 |0 |4 |4 |-8|4  |0  |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |-4|-4|0 |-4|0 |8 |4 |-4|0  |-8|4 |0 |-4|-4 |0  |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |4 |-4|8 |4 |0 |0 |4 |0 |-4 |4 |8 |-4|0 |0  |-4 |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |-4|0 |4 |4 |-8|4 |0 |-4|0  |4 |0 |8 |4 |0  |4  |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |0 |0 |0 |0 |0 |0 |0 |-4|4  |4 |-4|4 |-4|-4 |-12|
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |0 |-4|-4|0 |0 |-4|-4|-8|0  |-4|4 |0 |8 |4  |-4 |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |8 |-4|4 |-8|0 |4 |-4|4 |4  |0 |0 |4 |4 |0  |0  |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |8 |0 |-8|8 |0 |8 |0 |0 |0  |0 |0 |0 |0 |0  |0  |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |-4|8 |-4|-4|0 |4 |0 |4 |0  |4 |8 |0 |4 |0  |-4 |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |4 |4 |0 |-4|8 |0 |4 |-8|-4 |4 |0 |4 |0 |0  |4  |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |4 |4 |0 |-4|-8|0 |4 |-4|0  |0 |-4|-8|4 |-4 |0  |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+
|0 |-4|-8|-4|-4|0 |4 |0 |0 |-4 |8 |-4|-4|0 |4  |0  |
+--+--+--+--+--+--+--+--+--+---+--+--+--+--+---+---+

Walsh Spectrum representation (except first row and column):

.. image:: /images/Nibble.png
   :width: 750 px
   :align: center

`Linear Profile <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSub.lp>`_

`Differential Profile <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSub.dp>`_

`Autocorrelation Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSub.ac>`_

Other useful information in cryptanalysis
-----------------------------------------

Cycle structure:

+--------------+------------------+
| Cycle length | Number of cycles |
+==============+==================+
| 2            | 1                |
+--------------+------------------+
| 14           | 1                |
+--------------+------------------+

There are 7 linear structures:

.. code-block:: console

	([0 0 1 1],[0 1 1 1])
	([0 1 0 0],[0 1 0 1])
	([1 0 1 0],[1 1 1 1])
	([1 0 1 1],[0 1 0 1])
	([1 1 0 1],[1 0 0 1])
	([1 1 1 0],[1 1 1 0])
	([1 1 1 1],[0 1 0 1])

It has no fixed points and 2 negated fixed points: (0,0,1,0), (0,1,1,1)

NibbleSubInv
============

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSubInv.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSubInv.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSubInv.anf>`_

`Characteristic Function <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSubInv.char>`_

`Walsh Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSubInv.wal>`_

+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|16|0 |0  |0  |0 |0 |0 |0 |0  |0  |0 |0 |0 |0 |0  |0  |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |0 |0  |0  |4 |-4|4 |-4|0  |0  |8 |8 |-4|4 |4  |-4 |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |-4|-4 |0  |0 |-4|-4|0 |0  |-4 |-4|0 |8 |4 |4  |-8 |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |-4|-4 |0  |-4|0 |8 |4 |0  |-4 |4 |-8|-4|0 |0  |-4 |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |0 |0  |0  |-4|-4|4 |4 |0  |0  |-8|8 |-4|-4|-4 |-4 |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |0 |0  |0  |-8|0 |0 |-8|0  |0  |0 |0 |0 |8 |-8 |0  |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |-4|-4 |0  |-4|8 |0 |4 |0  |-4 |4 |8 |4 |0 |0  |4  |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |12|-4 |0  |0 |4 |4 |0 |0  |-4 |-4|0 |0 |4 |4  |0  |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |4 |0  |4  |0 |-4|0 |-4|-4 |-8 |4 |0 |4 |-8|-4 |0  |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |4 |0  |-12|-4|0 |-4|0 |4  |0  |4 |0 |0 |-4|0  |-4 |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |0 |4  |-4 |0 |-8|4 |4 |4  |-4 |0 |0 |4 |4 |0  |8  |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |0 |4  |-4 |4 |4 |8 |0 |-4 |4  |0 |0 |8 |0 |-4 |-4 |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |4 |0  |4  |4 |0 |-4|8 |4  |0  |4 |0 |0 |4 |-8 |-4 |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |4 |0  |4  |-8|-4|0 |4 |-4 |8  |4 |0 |4 |0 |4  |0  |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |0 |-12|-4 |4 |-4|0 |0 |-4 |4  |0 |0 |0 |0 |-4 |4  |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+
|0 |0 |4  |-4 |0 |0 |-4|4 |-12|-4 |0 |0 |-4 |4 |0 |0  |
+--+--+---+---+--+--+--+--+---+---+--+--+--+--+---+---+

Walsh Spectrum representation (except first row and column):

.. image:: /images/NibbleSubInv.png
   :width: 750 px
   :align: center

`Linear Profile <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSubInv.lp>`_

`Differential Profile <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSubInv.dp>`_

`Autocorrelation Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/NibbleSubInv.ac>`_

Other useful information in cryptanalysis
-----------------------------------------

Cycle structure:

+--------------+------------------+
| Cycle length | Number of cycles |
+==============+==================+
| 2            | 1                |
+--------------+------------------+
| 14           | 1                |
+--------------+------------------+

There are 7 linear structures:

.. code-block:: console

	([0 0 1 0],[0 0 1 0])
	([0 1 0 1],[1 0 1 1])
	([1 0 0 0],[0 0 1 1])
	([1 0 0 0],[1 0 0 0])
	([1 0 0 0],[1 0 1 1])
	([1 0 1 0],[0 0 0 1])
	([1 1 0 1],[1 0 1 1])

It has no fixed points and 2 negated fixed points: (1,0,0,0), (1,1,0,1)

MixColumn
=========

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/mixcolumn.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/mixcolumn.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/mixcolumn.anf>`_

`Characteristic Function <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/mixcolumn.char>`_

`Walsh Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/mixcolumn.wal>`_

Walsh Spectrum representation (except first row and column):

.. image:: /images/mixcolumn.png
   :width: 750 px
   :align: center

`Linear Profile <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/mixcolumn.lp>`_

`Differential Profile <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/mixcolumn.dp>`_

`Autocorrelation Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/mixcolumn.ac>`_

Other useful information in cryptanalysis
-----------------------------------------

Cycle structure:

+--------------+------------------+
| Cycle length | Number of cycles |
+==============+==================+
| 1            | 16               |
+--------------+------------------+
| 2            | 120              |
+--------------+------------------+

There 65025 linear structures

It has 15 fixed points: (0,0,0,0,0,0,0,0), (0,0,0,1,0,0,0,1), (0,0,1,0,0,0,1,0), (0,0,1,1,0,0,1,1), (0,1,0,0,0,1,0,0), (0,1,0,1,0,1,0,1), (0,1,1,0,0,1,1,0), (0,1,1,1,0,1,1,1), (1,0,0,0,1,0,0,0), (1,0,0,1,1,0,0,1), (1,0,1,0,1,0,1,0), (1,0,1,1,1,0,1,1), (1,1,0,0,1,1,0,0), (1,1,0,1,1,1,0,1), (1,1,1,0,1,1,1,0)

It has 16 negated fixed points: (0,0,0,0,1,1,1,0), (0,0,0,1,1,1,1,1), (0,0,1,0,1,1,0,0), (0,0,1,1,1,1,0,1), (0,1,0,0,1,0,1,0), (0,1,0,1,1,0,1,1), (0,1,1,0,1,0,0,0), (0,1,1,1,1,0,0,1), (1,0,0,0,0,1,1,0), (1,0,0,1,0,1,1,1), (1,0,1,0,0,1,0,0), (1,0,1,1,0,1,0,1), (1,1,0,0,0,0,1,0), (1,1,0,1,0,0,1,1), (1,1,1,0,0,0,0,0), (1,1,1,1,0,0,0,1)

ks0
===

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/ks0.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/ks0.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/ks0.anf>`_

`Walsh Spectrum (each row represents a column of Walsh Spectrum) <https://github.com/jacubero/VBF/blob/master/miniAES/ks0.wal.gz>`_

`Linear Profile (each row represents a column of Linear Profile) <https://github.com/jacubero/VBF/blob/master/miniAES/ks0.lp.gz>`_

Other useful information in cryptanalysis
-----------------------------------------

Cycle structure:

+--------------+------------------+
| Cycle length | Number of cycles |
+==============+==================+
| 1            | 65536            |
+--------------+------------------+

ks1
===

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/ks1.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/ks1.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/ks1.anf>`_

`Walsh Spectrum (each row represents a column of Walsh Spectrum) <https://github.com/jacubero/VBF/blob/master/miniAES/ks1.wal.gz>`_

`Linear Profile (each row represents a column of Linear Profile) <https://github.com/jacubero/VBF/blob/master/miniAES/ks1.lp.gz>`_

Other useful information in cryptanalysis
-----------------------------------------

Cycle structure:

+--------------+------------------+
| Cycle length | Number of cycles |
+==============+==================+
| 1            | 1                |
+--------------+------------------+
| 5            | 2                |
+--------------+------------------+
| 10           | 2                |
+--------------+------------------+
| 28           | 1                |
+--------------+------------------+
| 60           | 1                |
+--------------+------------------+
| 1223         | 1                |
+--------------+------------------+
| 26097        | 1                |
+--------------+------------------+
| 38097        | 1                |
+--------------+------------------+

ks2
===

Representations
---------------

Polynomial representation in ANF:

`f1 <https://github.com/jacubero/VBF/blob/master/miniAES/f1.pdf>`_

`f2 <https://github.com/jacubero/VBF/blob/master/miniAES/f2.pdf>`_

`f3 <https://github.com/jacubero/VBF/blob/master/miniAES/f3.pdf>`_

`f4 <https://github.com/jacubero/VBF/blob/master/miniAES/f4.pdf>`_

`f5 <https://github.com/jacubero/VBF/blob/master/miniAES/f5.pdf>`_

`f6 <https://github.com/jacubero/VBF/blob/master/miniAES/f6.pdf>`_

`f7 <https://github.com/jacubero/VBF/blob/master/miniAES/f7.pdf>`_

`f8 <https://github.com/jacubero/VBF/blob/master/miniAES/f8.pdf>`_

`f9 <https://github.com/jacubero/VBF/blob/master/miniAES/f9.pdf>`_

`f10 <https://github.com/jacubero/VBF/blob/master/miniAES/f10.pdf>`_

`f11 <https://github.com/jacubero/VBF/blob/master/miniAES/f11.pdf>`_

`f12 <https://github.com/jacubero/VBF/blob/master/miniAES/f12.pdf>`_

`f13 <https://github.com/jacubero/VBF/blob/master/miniAES/f13.pdf>`_

`f14 <https://github.com/jacubero/VBF/blob/master/miniAES/f14.pdf>`_

`f15 <https://github.com/jacubero/VBF/blob/master/miniAES/f15.pdf>`_

`f16 <https://github.com/jacubero/VBF/blob/master/miniAES/f16.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/ks2.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/ks2.anf>`_

`Walsh Spectrum (each row represents a column of Walsh Spectrum) <https://github.com/jacubero/VBF/blob/master/miniAES/ks2.wal.gz>`_

`Linear Profile (each row represents a column of Linear Profile) <https://github.com/jacubero/VBF/blob/master/miniAES/ks2.lp.gz>`_

Other useful information in cryptanalysis
-----------------------------------------

Cycle structure:

+--------------+------------------+
| Cycle length | Number of cycles |
+==============+==================+
| 1            | 1                |
+--------------+------------------+
| 12           | 1                |
+--------------+------------------+
| 15           | 3                |
+--------------+------------------+
| 30           | 1                |
+--------------+------------------+
| 109          | 1                |
+--------------+------------------+
| 385          | 1                |
+--------------+------------------+
| 831          | 1                |
+--------------+------------------+
| 2472         | 1                |
+--------------+------------------+
| 3617         | 1                |
+--------------+------------------+
| 9775         | 1                |
+--------------+------------------+
| 16777        | 1                |
+--------------+------------------+
| 31482        | 1                |
+--------------+------------------+

mini-AES
========

Algebraic degree from key 00000 to 65535 is equal to 14

`Cycle structure from key 00000 to 65535 <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/cycle.pdf>`_

`Fixed and negated points from key 00000 to 65535 <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/points.pdf>`_

`Nonlinearities from key 00000 to 65535 <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/nl.pdf>`_

`Nonlinearities in ascendent order <https://raw.githubusercontent.com/jacubero/VBF/master/miniAES/fi.s>`_

Graphical display of the distribution of the nonlinearities of mini-AES:

.. image:: /images/hist-miniAES.jpeg
   :width: 750 px
   :align: center

+---------------------------------------------------+
| Descriptive Statistics of mini-AES nonlinearities |
+====================+==============================+
| Unique Values      | 130                          |
+--------------------+------------------------------+
| Min                | 31432                        |
+--------------------+------------------------------+
| Max                | 32040                        |
+--------------------+------------------------------+
| Mean               | 31912.9894                   |
+--------------------+------------------------------+
| Mean Deviation     | 8.6571                       |
+--------------------+------------------------------+
| 1st Quartile       | 31880                        |
+--------------------+------------------------------+
| Median             | 31924                        |
+--------------------+------------------------------+
| 3rd Quartile       | 31960                        |
+--------------------+------------------------------+
| Mode               | 31952                        |
+--------------------+------------------------------+
| Range              | 608                          |
+--------------------+------------------------------+
| Variance           | 3903.8642                    |
+--------------------+------------------------------+
| Standard Deviation | 62.4809                      |
+--------------------+------------------------------+
| Kkewness           | -1.092059                    |
+--------------------+------------------------------+
| Kurtosis           | 1.79284                      |
+--------------------+------------------------------+
| P0.5               | 31692                        |
+--------------------+------------------------------+
| P1                 | 31720                        |
+--------------------+------------------------------+
| P5                 | 31796                        |
+--------------------+------------------------------+
| P95                | 31992                        |
+--------------------+------------------------------+
| P99                | 32012                        |
+--------------------+------------------------------+
| P99.5              | 32016                        |
+--------------------+------------------------------+

