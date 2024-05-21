from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Python package for calculating three TME-based scores'
LONG_DESCRIPTION = 'This is the Python package for the ESTIMATE algorithm, ISTMEscore and SIA score. The first two modules are the adapted versions of the R package `estimate` and `ISTMEscore`. Link for papers: '

# Setting up
setup(
        name="TMEscore", 
        version=VERSION,
        author="Qilu Zhou",
        author_email="<qiluzhou@umass.edu>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'TME score'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
           # "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)