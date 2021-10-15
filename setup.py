from setuptools import setup

setup(name='terra',
      version='0.2',
      description='Script to determine missing rocky material of a star',
      url='https://github.com/ramstojh/terra',
      author='Jhon Yana',
      author_email='ramstojh@alumni.usp.br',
      license='MIT',
      packages=['terra'],
      install_requires=['numpy', 'pandas', 'tqdm', 'astropy', 'matplotlib'],
      zip_safe=False)
