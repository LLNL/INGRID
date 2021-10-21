*************************
The INGRID parameter file
*************************

Background
==========

INGRID has been designed to be controlled from a single configuration/parameter file when operating in GUI mode. This specially formatted `YAML <https://yaml.org/spec/1.2/spec.html#Introduction>`_ file is similar to Fortran namelist files due to the `key`-`value` structure it contains. 

Here is a snippet of what a YAML formatted file can look like.

.. code-block:: YAML

    # Comments are supported and follow python convention!
    # YAML entries are a mapping of the form:
    #
    #   key: value
    #
    # It follows that python interprets the YAML file as a dict

    my_YAML_entry:
        # YAML file use spaces as indentation.
        # 2 or 4 spaces (pick one and stick to it) indicate a nested item.
        # Here we use 4 spaces.
        my_str_key: 'Hello, YAML!'  # my_str_key = 'Hello, YAML!'
        my_int_key: 32  # my_int_key = 32
        some_float: 3.14159  # some_float = 3.14159
        my_true_bool: true # bool case-insensitive 
        my_false_bool: 0  # '1' '0' bool values supported.

        # Empty lines within file are ok.

        # Just remember it's the spaces that matter.

        my_sub_dict: # Nested dicts are supported.
            another_key: ending my example code-block here. 
            # The above value entry will be interpreted as a string
            # (no quotations needed)

Although the INGRID parameter file contains many controls for a user, it does not stray from the patterns illustrated above (no advanced YAML knowledge required).

.. tip:: While operating INGRID in GUI mode, keep your favorite text-editor handy with the parameter-file in use loaded. You will be making frequent edits to this core parameter/configuration file.

Although not necessary for following tutorial, detailed documentation of the INGRID parameter-file can be found `here <parameter_file_documentation>`_.


A single-null example file
==========================

Users of INGRID have plenty of controls available for fine-tuning their final grid. This section explains how to navigate some key controls within the INGRID parameter file.

We will be using the pre-populated INGRID parameter file ``DIIID_SNL.yml`` from the ``example_files`` directory in the single-null walk-through. Open this file in your preferred text-editor. At the time of writing this documentation, the parameter file ``DIIID_SNL.yml`` contains the following:

.. code-block:: YAML

    # ---------------------------------------------------
    # User data directories
    # ---------------------------------------------------
    dir_settings:
      eqdsk: ../data/SNL/DIII-D/  # dir containing eqdsk
      limiter: .  # dir containing limiter
      patch_data: ../data/SNL/DIII-D/  # dir containing patch data
      target_plates: ../data/SNL/DIII-D/ # dir containing target plates

    # ---------------------------------------------------
    # eqdsk file name
    # ---------------------------------------------------
    eqdsk: neqdsk

    # ---------------------------------------------------
    # General grid settings
    # ---------------------------------------------------
    grid_settings:
      # ----------------------------------------------------------------------------
      # Settings for grid generation (num cells, transforms, skewness_correction)
      # ----------------------------------------------------------------------------
      grid_generation:
        skewness_correction:
          all:
            active: True # true, 1 also valid.
            resolution: 1000
            theta_max: 120.0
            theta_min: 80.0
        np_default: 3
        nr_default: 3
        poloidal_f_default: x, x
        radial_f_default: x, x
      # ---------------------------------------------------
      # guard cell size
      # ---------------------------------------------------
      guard_cell_eps: 0.001
      # ---------------------------------------------------
      # num levels in efit plot
      # ---------------------------------------------------
      nlevs: 30
      # ---------------------------------------------------
      # num xpts
      # ---------------------------------------------------
      num_xpt: 1
      patch_generation:
        strike_pt_loc: target_plates # 'limiter' or 'target_plates'
        rmagx_shift: 0.0
        zmagx_shift: 0.0
      # ---------------------------------------------------
      # Psi levels
      # ---------------------------------------------------
      psi_1: 1.066
      psi_core: 0.95
      psi_pf_1: 0.975
      # ---------------------------------------------------
      # magx coordinates
      # ---------------------------------------------------
      rmagx: 1.75785604
      zmagx: -0.0292478683
      # ---------------------------------------------------
      # xpt coordinates
      # ---------------------------------------------------
      rxpt: 1.300094032687
      zxpt: -1.133159375302
      # ---------------------------------------------------
      # Filled contours vs contour lines
      # ---------------------------------------------------
      view_mode: filled

    # ---------------------------------------------------
    # Saved patch settings
    # ---------------------------------------------------
    patch_data:
      file: LSN_patches_1597099640.npy
      preferences:
        new_file: true
        new_fname: LSN_patches_1597099640.npy
      use_file: false

    # ---------------------------------------------------
    # Integrator
    # ---------------------------------------------------
    integrator_settings:
      dt: 0.01
      eps: 5.0e-06
      first_step: 5.0e-05
      max_step: 0.064
      step_ratio: 0.02
      tol: 0.005

    # ---------------------------------------------------
    # Limiter settings
    # ---------------------------------------------------
    limiter:
      file: ''
      use_efit_bounds: false

    # ---------------------------------------------------
    # target plate settings
    # ---------------------------------------------------
    target_plates:
      plate_E1:
        file: d3d_otp.txt
        zshift: -1.6
      plate_W1:
        file: d3d_itp.txt
        zshift: -1.6

Let's highlight some important entries that are often used when operating INGRID for single-null cases (basic usage). Advanced tutorials will also be provided.

Setting data paths
==================

A user can provide a string that indicates the path to certain data. This is used to tell INGRID where to look for EFIT data, target plate coordinate, limiter coordinates, and patch-data (for reconstruction). We can set these paths by editing the entry:

.. code-block:: YAML

    # ---------------------------------------------------
    # User data directories
    # ---------------------------------------------------
    dir_settings:
      eqdsk: ../data/SNL/DIII-D/  # dir containing eqdsk
      limiter: .  # dir containing limiter
      patch_data: ../data/SNL/DIII-D/  # dir containing patch data
      target_plates: ../data/SNL/DIII-D/ # dir containing target plates

.. note:: INGRID supports both absolute paths and paths relative to where INGRID has been launched.

If ``dir_settings`` is missing any entries, INGRID will (internally) set the missing values to a default value of ``'.'`` (current working directory). This holds even if ``dir_settings`` is omitted from the parameter file.

.. note:: ``dir_settings`` entries are the **directory** to look for data and NOT the file itself.

Providing an EQDSK file
=======================

The user provides the actual EQDSK file name separate from the ``dir_settings`` entry. We provide this at the global YAML level under entry ``eqdsk``. That is:

.. code-block:: YAML

    # ---------------------------------------------------
    # eqdsk file name
    # ---------------------------------------------------
    eqdsk: neqdsk

.. note:: In this example, INGRID searches for the file ``neqdsk`` within the directory ``../data/SNL/DIII-D/`` (relative to the launch point) since ``dir_settings['eqdsk']`` was set to ``../data/SNL/DIII-D/`` (see above).

Defining target plates
=======================

All target plate settings are under the global INGRID parameter file entry ``target_plates``. We see this as:

.. code-block:: YAML

    # ---------------------------------------------------
    # target plate settings
    # ---------------------------------------------------
    target_plates:
      plate_E1:
        file: d3d_otp.txt
        zshift: -1.6
      plate_W1:
        file: d3d_itp.txt
        zshift: -1.6

INGRID adopts a N-S-E-W compass direction notation in order to help generalize and simplify grid generation. It is important for a user to eventually learn these conventions. A detailed discussion of INGRID's naming conventions can be found `here <ingrid_notation>`__.

For now (in the case of a lower single-null configuration), note that entries ``plate_E1`` and ``plate_W1`` correspond to the `outer` and `inner` target plates, respectively. Each plate entry recognizes sub-entries ``file`` (file name to load), ``zshift`` (z-translation) and ``rshift`` (r-translation, not utilized and internal to INGRID defaults to ``0.0``). 


Defining x-points, magnetic-axis, and psi-levels
=================================================

Settings for x-point coordinates, magnetic-axis coordinates, and psi-levels are found under the global INGRID parameter file entry ``grid_settings``.

.. code-block:: YAML

    # ---------------------------------------------------
    # General grid settings
    # ---------------------------------------------------
    grid_settings:
      # ...
      # ...
      # Other items currently not of interest...
      # ...
      # ...

      # ---------------------------------------------------
      # num xpts
      # ---------------------------------------------------
      num_xpt: 1

      # ---------------------------------------------------
      # Psi levels
      # ---------------------------------------------------
      psi_1: 1.066  # SOL
      psi_core: 0.95  # CORE
      psi_pf_1: 0.975  # PRIVATE-FLUX

      # ---------------------------------------------------
      # magx coordinates
      # ---------------------------------------------------
      rmagx: 1.75785604
      zmagx: -0.0292478683

      # ---------------------------------------------------
      # primary xpt coordinates
      # ---------------------------------------------------
      rxpt: 1.300094032687
      zxpt: -1.133159375302

.. warning:: The entry ``num_xpt`` is one of the most important entries in the INGRID parameter file since it determines INGRID's method of analysis. Dealing with more than one x-point requires a more in-depth understanding of the parameter file, so ensure this is set to the correct number of x-points.

