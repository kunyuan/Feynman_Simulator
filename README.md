# Feynman Simulator
A generic purpose quantum simulator using Feynman diagram technique

## Import project to CLion
   To create a CLion project, you have to do the followings:
   * _File | Import Project_, and choose __src__ folder, then in the dialog window, choose "Open Project"
     instead of "Rewrite CMakeList".
   * The last step will make __src__ folder as the project root (By default, CLion expect a CMakeLists.txt file to be provided in the root
     directory of all project-related files). Use _Tools | CMake | Change Project Root_ to change the root
     to the entire codebase.
   * In _Run | Edit Configurations_, change the Target "simulator.exe"'s _Program arguments_ to _-p 0_
     (or other input parameters), and _Working directory_ to the foler contains "simulator.exe" file. 
   * Click _Run | build_, you should be able to see the "simulator.exe" in the code base folder.

## User tips

1. Simulation status check

   To access remote notebook on server
   run __ipython notebook --no-browser --ip=* --port=8889__,
   then use: __http://...:8889__ to access notebook.
