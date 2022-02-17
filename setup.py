from setuptools import setup, Extension

setup (name = 'spkmeansmodule',
       version = '1.0',
       description = 'allows the use of wanted *.c files',
       ext_modules =[
           Extension(
               'spkmeansmodule',
                sources=['spkmeans.c','spkmeansmodule.c']
           )
       ]
)