from setuptools import setup, find_packages
 
classifiers = [
  'Development Status :: 5 - Production/Stable',
  'Intended Audience :: Education',
  'Operating System :: Microsoft :: Windows :: Windows 10',
  'License :: OSI Approved :: MIT License',
  'Programming Language :: Python :: 3'
]
 
setup(
  name='FrozenYoghourt',
  version='0.0.5',
  description='Numerical and Symbolic Manipulation for Quantum Computing',
  long_description=open('README.txt').read() + '\n\n' + open('CHANGELOG.txt').read(),
  url='',  
  author='Peter Montgomery',
  author_email='petermontgomery056@gmail.com',
  license='MIT', 
  classifiers=classifiers,
  keywords='quantum computing', 
  packages=find_packages(),
  install_requires=[''] 
)