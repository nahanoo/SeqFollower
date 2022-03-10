
from setuptools import setup, find_packages

setup(name='insertion_detection',
      version='1.0',
      description="Detects and visualizes large insertions in mutated bacterial strains.",
      author='Eric Ulrich',
      url='https://github.com/nahanoo/insertion_detection',
      packages=['insertion_detection'],
      install_requires=['pandas',
                        'pysam',
                        'Bio',
                        'dna_features_viewer'],
      entry_points={
          'console_scripts': [
              'detect_insertions = insertion_detection.main:main'
          ]
      }
      )
