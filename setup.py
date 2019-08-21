from setuptools import setup, find_packages, Extension

setup(name='aquila_umap',
      version='1.0',
      description='this is a program to generate Umap for Aquila diploid assembly',
      author='xzhou',
      author_email='xzhou15@cs.stanford.edu',
      packages=['aquila_umap',],
      entry_points={'console_scripts':['aquila_umap=aquila_umap.Aquila_Umap:main']},
      zip_safe=False)
