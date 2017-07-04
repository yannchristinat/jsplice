from setuptools import setup

setup(name='jsplice',
      version='1.0.3',
      description='jSplice - Rapid detection of differential alternative splicing events from RNA-seq',
      url='https://github.com/yannchristinat/jsplice',
      author='Yann Christinat',
      author_email='yann.christinat@gmail.com',
      license='GNU GPL 3.0',
      packages=['jsplice','jsplice.lib'],
	  install_requires=['scipy','numpy'],
	  scripts=['bin/jsplice','bin/starJxn2bed'],
      zip_safe=False)