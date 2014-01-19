viennafem-dev
==============

Developer repository for ViennaFEM. Visit http://viennamath.sourceforge.net/ if you are looking for the latest releases.

How to build:

 - Clone this repository.
 - Clone the developer repositories of [ViennaCL](https://github.com/viennacl/viennacl-dev/), [ViennaData](https://github.com/viennadata/viennadata-dev/), [ViennaGrid](https://github.com/viennagrid/viennagrid-dev/), and [ViennaMath](https://github.com/viennamath/viennamath-dev/) into separate folders. If you know how to do this automatically with CMake, please let us know.
 - Copy the source folders viennacl/, viennadata/, viennagrid/, and viennamath/ from the respective project to the main ViennaFEM folder (i.e. where the folder viennafem/ resides). Feel free to create symbolic links instead of copying.
 - Create the build directory and build using CMake:

        $> mkdir build/
        $> cd build/
        $ build> cmake ..
        $ build> make

 - Run the examples from the build folder, e.g.

        $ build> examples/tutorials/poisson_2d



