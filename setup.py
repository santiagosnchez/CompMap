import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='CompMap',  
     version='1.0',
     scripts=['CompMap.py'] ,
     author="Santiago Sanchez-Ramirez",
     author_email="santiago.snchez@gmail.com",
     description="Competitive read-mapping for allele-specific expression read-counting",
     long_description=long_description,
   long_description_content_type="text/markdown",
     url="https://github.com/santiagosnchez/CompMap",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: GPL3 License",
         "Operating System :: OS Independent",
     ],
 )
