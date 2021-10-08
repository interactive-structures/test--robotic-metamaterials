# Robotic Metamaterials
 
 TODO: update and integrate contents from previous README

Additional instructions: I use Intel MKL, instruct MSVC in VS project properties to use it (needs to be installed and added to Path, TODO: copy contents of path for reference!)

---------------------------------

git submodule to specific commit

1) add submodule:
```
git submodule add <URL> <path>
```
example: 
` git submodule add https://github.com/libigl/libigl.git ./dependencies/libigl `

2) checkout the tag/commit we should track. note it is not encoded explicitly but apparently is tracked by the main project's git config. havent tested it

```
git checkout tags/<tagname> -b <tagname>-branch
```
example: 
` git checkout tags/v2.3.0 -b v2.3.0-branch `

this switches to the tag commit and creates a new branch of that commit