***
DES
***

.. _secExamplesDES

Description
===========

*DES* is a block cipher used for securing sensitive but unclassified material by U.S. Government agencies and became the de facto encryption standard for commercial transactions in the private sector up to December 6, 2001. It was superseded by AES. DES has eight 6x4 S-boxes: S1, S2, S3, S4, S5, S6, S7, S8.

Summary
=======

+-------+------+-------+------+-------+------+---------+----------------+------------+------+
| S-box | *NL* | *NL2* | *LD* | *DEG* | *AI* | *MAXAC* | :math:`\sigma` | *LP*       | *DP* |
+=======+======+=======+======+=======+======+=========+================+============+======+
| S1    | 14   | 8     | 4    | 4     | 2    | 48      | 36736          | 0.31640625 | 0.25 |
+-------+------+-------+------+-------+------+---------+----------------+------------+------+
| S2    | 16   | 10    | 2    | 4     | 3    | 56      | 25984          | 0.25       | 0.25 |
+-------+------+-------+------+-------+------+---------+----------------+------------+------+
| S3    | 16   | 10    | 4    | 4     | 3    | 48      | 24064          | 0.25       | 0.25 |
+-------+------+-------+------+-------+------+---------+----------------+------------+------+
| S4    | 16   | 8     | 0    | 3     | 2    | 64      | 40960          | 0.25       | 0.25 |
+-------+------+-------+------+-------+------+---------+----------------+------------+------+
| S5    | 12   | 10    | 6    | 4     | 2    | 40      | 47104          | 0.390625   | 0.25 |
+-------+------+-------+------+-------+------+---------+----------------+------------+------+
| S6    | 18   | 10    | 4    | 5     | 3    | 48      | 19456          | 0.19140625 | 0.25 |
+-------+------+-------+------+-------+------+---------+----------------+------------+------+
| S7    | 14   | 10    | 4    | 5     | 3    | 48      | 34048          | 0.31640625 | 0.25 |
+-------+------+-------+------+-------+------+---------+----------------+------------+------+
| S8    | 16   | 10    | 4    | 4     | 3    | 48      | 24064          | 0.25       | 0.25 |
+-------+------+-------+------+-------+------+---------+----------------+------------+------+

S1
==

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S1/S1.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S1/S1.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S1/S1.anf>`_

`Characteristic function <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S1/S1.char>`_

`Walsh Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S1/S1.wal>`_

Walsh Spectrum representation (except first row and column):

.. image:: /images/S1WS.png
   :width: 750 px
   :align: center

`Linear Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S1/S1.lp>`_

`Differential Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S1/S1.dp>`_

`Autocorrelation Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S1/S1.ac>`_

Other useful information in cryptanalysis
-----------------------------------------

There are no linear structures

S2
==

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S2/S2.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S2/S2.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S2/S2.anf>`_

`Characteristic function <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S2/S2.char>`_

`Walsh Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S2/S2.wal>`_

Walsh Spectrum representation (except first row and column):

.. image:: /images/S2WS.png
   :width: 750 px
   :align: center

`Linear Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S2/S2.lp>`_

`Differential Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S2/S2.dp>`_

`Autocorrelation Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S2/S2.ac>`_

Other useful information in cryptanalysis
-----------------------------------------

There are no linear structures

S3
==

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S3/S3.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S3/S3.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S3/S3.anf>`_

`Characteristic function <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S3/S3.char>`_

`Walsh Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S3/S3.wal>`_

Walsh Spectrum representation (except first row and column):

.. image:: /images/S3WS.png
   :width: 750 px
   :align: center

`Linear Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S3/S3.lp>`_

`Differential Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S3/S3.dp>`_

`Autocorrelation Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S3/S3.ac>`_

Other useful information in cryptanalysis
-----------------------------------------

There are no linear structures

S4
==

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S4/S4.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S4/S4.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S4/S4.anf>`_

`Characteristic function <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S4/S4.char>`_

`Walsh Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S4/S4.wal>`_

Walsh Spectrum representation (except first row and column):

.. image:: /images/S4WS.png
   :width: 750 px
   :align: center

`Linear Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S4/S4.lp>`_

`Differential Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S4/S4.dp>`_

`Autocorrelation Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S4/S4.ac>`_

Other useful information in cryptanalysis
-----------------------------------------

There are 9 linear structures:

.. code-block:: console

   ([0 0 0 0 0 1],[0 0 1 1])
   ([0 0 0 0 0 1],[1 1 0 0])
   ([0 0 0 0 0 1],[1 1 1 1])
   ([1 0 1 1 1 0],[0 1 0 1])
   ([1 0 1 1 1 0],[1 0 1 0])
   ([1 0 1 1 1 0],[1 1 1 1])
   ([1 0 1 1 1 1],[0 1 1 0])
   ([1 0 1 1 1 1],[1 0 0 1])
   ([1 0 1 1 1 1],[1 1 1 1])

S5
==

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S5/S5.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S5/S5.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S5/S5.anf>`_

`Characteristic function <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S5/S5.char>`_

`Walsh Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S5/S5.wal>`_

Walsh Spectrum representation (except first row and column):

.. image:: /images/S5WS.png
   :width: 750 px
   :align: center

`Linear Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S5/S5.lp>`_

`Differential Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S5/S5.dp>`_

`Autocorrelation Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S5/S5.ac>`_

Other useful information in cryptanalysis
-----------------------------------------

There are no linear structures

S6
==

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S6/S6.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S6/S6.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S6/S6.anf>`_

`Characteristic function <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S6/S6.char>`_

`Walsh Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S6/S6.wal>`_

Walsh Spectrum representation (except first row and column):

.. image:: /images/S6WS.png
   :width: 750 px
   :align: center

`Linear Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S6/S6.lp>`_

`Differential Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S6/S6.dp>`_

`Autocorrelation Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S6/S6.ac>`_

Other useful information in cryptanalysis
-----------------------------------------

There are no linear structures

S7
==

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S7/S7.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S7/S7.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S7/S7.anf>`_

`Characteristic function <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S7/S7.char>`_

`Walsh Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S7/S7.wal>`_

Walsh Spectrum representation (except first row and column):

.. image:: /images/DES-S7WS.png
   :width: 750 px
   :align: center

`Linear Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S7/S7.lp>`_

`Differential Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S7/S7.dp>`_

`Autocorrelation Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S7/S7.ac>`_

Other useful information in cryptanalysis
-----------------------------------------

There are no linear structures

S8
==

Representations
---------------

`Polynomial representation in ANF <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S8/S8.pdf>`_

`Truth Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S8/S8.tt>`_

`ANF Table <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S8/S8.anf>`_

`Characteristic function <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S8/S8.char>`_

`Walsh Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S8/S8.wal>`_

Walsh Spectrum representation (except first row and column):

.. image:: /images/S8WS.png
   :width: 750 px
   :align: center

`Linear Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S8/S8.lp>`_

`Differential Profile <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S8/S8.dp>`_

`Autocorrelation Spectrum <https://raw.githubusercontent.com/jacubero/VBF/master/DES/S8/S8.ac>`_

Other useful information in cryptanalysis
-----------------------------------------

There are no linear structures
