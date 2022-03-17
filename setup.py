from setuptools import setup

setup(name='SeqFollower',
      version='1.0',
      description='Detect inserted or deleted sequences from assemblies of evolved bacterial strains.',
      author='Eric Ulrich',
      url='https://github.com/nahanoo/SeqFollower',
      packages=['deletion_detection','insertion_detection','hgt_detection'],
      install_requires=['pandas',
                        'pysam',
                        'Bio',
                        'dna_features_viewer'],
      entry_points={
          'console_scripts': [
              'detect_deletions = deletion_detection.main:main',
              'detect_insertions = insertion_detection.main:main',
              'detect_hgts = hgt_detection.main:main',
          ]
      }
     )
