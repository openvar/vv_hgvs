from setuptools import setup, find_packages
from sys import version_info

short_description = "HGVS Parser, Formatter, Mapper, Validator"
with open("doc/description.txt") as f:
    long_description = f.read()

if version_info < (3, ):
    version_specific_requirements = ['unicodecsv']
else:
    version_specific_requirements = []


setup(license="Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0)",
      long_description=long_description,
      use_scm_version=True,
      zip_safe=True,
      author="HGVS Contributors",
      author_email = 'hgvs-discuss@gmail.com',
      description = short_description.replace("\n", " "),
      name="vvhgvs",
      package_data={"vvhgvs": ["_data/*"]},
      packages=find_packages(),
      url="https://github.com/openvar/vv_hgvs",
      entry_points={
          'console_scripts': [
              'vvhgvs-shell = vvhgvs.shell:shell'
          ],
      },
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "Intended Audience :: Developers",
          "Intended Audience :: Healthcare Industry",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: Apache Software License",
          "Operating System :: OS Independent",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 2",
          "Programming Language :: Python",
          "Topic :: Database :: Front-Ends",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Topic :: Scientific/Engineering :: Medical Science Apps.",
      ],
      keywords=[
          "bioinformatics",
          "computational biology",
          "genome variants",
          "genome variation",
          "genomic variants",
          "genomic variation",
          "genomics",
          "vvhgvs",
      ],
      install_requires=[
          "attrs>=17.4.0",  # https://github.com/biocommons/hgvs/issues/473
          "biocommons.seqrepo",
          "numpy",
          "biopython==1.74",    # 1.70 fails on rtd due to numpy absence
          "bioutils>=0.2.2",
          "configparser>=3.3.0",
          "enum34",
          "ipython<6",          # for vvhgvs-shell; >=6 for Py3 only
          "parsley",
          "psycopg2-binary",
          "six",
      ] + version_specific_requirements,
      setup_requires=[
          "pytest-runner",
          "setuptools_scm",
          "wheel",
      ],
      tests_require=[
          "pytest",
          "pytest-cov",
      ],
)


# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
