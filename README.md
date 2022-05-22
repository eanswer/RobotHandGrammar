# Hand Design System

#### Tested Platform

- Ubuntu 16.04, 18.04
- Windows 10 + Visual Studio 2019


#### How to run code?

###### Ubuntu
- Use CMake to compile the code
  ```
  cd RobotHandGrammar
  mkdir build
  cd build
  cmake ../ -DCMAKE_BUILD_TYPE=Release
  make
  ```
- Run the code
  ```
  cd build
  ./src/design
  ```
###### Windows 10

- *CMake*: Install latest version of cmake from [link](https://cmake.org/download/). Select the version 3.19.2, binary distribution  > Windows win64-x64 ZIP
- *g++*: check whether your compute has installed g++ by opening the terminal and run `g++ --version`. If it shows some version information, then you are done with this step, otherwise you can install g++ from [mingw](https://sourceforge.net/projects/mingw-w64).
- *Visual Studio 2019*: download and install visual studio 2019 from the [official website](https://visualstudio.microsoft.com/zh-hans/downloads/). The community version is free. During installation, make sure to select the *"Desktop development with C++"* option in the "Workloads" tab. If you forget to do so during installation, you can also do it after installation by referring the solution in [link](https://stackoverflow.com/questions/51668676/cmake-visual-studio-15-2017-could-not-find-any-instance-of-visual-studio).
- *clone the code", pick a folder you want to put the code in, open terminal in this folder, and run `git clone https://github.com/eanswer/RoboticHandDesign.git`. 
- *Configure the project*: 
  - open *cmake-gui* from the menu. 
  - Select the `RoboticHandDesign/Design/core` folder as *where is the source code*, and set  `RoboticHandDesign/Design/core/build` as *where to build the binaries* (create the folder if it is not there).
  - press *Configure*. Once it finishes, press *Generate*.
- *Run the code*:
  - Enter `RoboticHandDesign/Design/core/build`, double click `HandDesign.sln` to open the project in Visual Studio.
  - Find *Design* project in the *Solution Explorer* Panel, right click and select *Set as Startup Project".
  - Switch the mode from *Debug* to *Release* in the menu.
  - Now, if you press *Ctrl + F5*, you are able to run the code.

#### How does the UI work
- This UI provide the process to construct the palm structure by palm gramar, construct the finger structure by finger grammar, and change the morphology of the constructed palm and finger by cage-based deformation [[Xu et al. 2021]](http://diffhand.csail.mit.edu/)
- The work flow is: Palm Grammar -> Finger Grammar -> Cage-based Deformation
- Once you press the finish in one stage, the system will enter the next stage.
- During the Palm Grammar and Finger Grammar stages, the applicable rules are displayed on the left panel.
- In the Cage-based Deformation, click the component on the object while pressing Ctrl button to select the component to change the shape. If the selected shape is allowed change its shape, its changeable parameters (i.e. cage parameters) are displayed on the left panel and you are able to drag there to change the values. If there is nothing shown on the left panel, it means the geometry of the selected part is not modifiable (either because of manufacture constraints or constrained by the shapes of other parts).

#### Citation

If you find our paper or the code is useful, please consider citing:
```
@inproceedings{Zlokapa2022hand,
  title={An Integrated Design Pipeline for Tactile Sensing Robotic Manipulators},
  author={Zlokapa, Lara and Luo, Yiyue and Xu, Jie and Foshey, Michael and Wu, Kui and Agrawal Pulkit and Matusik, Wojciech},
  booktitle={2022 International conference on robotics and automation (ICRA)},
  year={2022},
  organization={IEEE}
}

```
