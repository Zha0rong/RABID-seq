from setuptools import setup, find_packages

try:
    from pip._internal.req import parse_requirements
except ImportError:
    from pip.req import parse_requirements

reqs = parse_requirements('requirements.txt', session=False)
requirements = [str(ir.req) for ir in reqs]

setup(
    name="rabid_seq",
    version="1.0.0",
    packages=find_packages(),
    author='Andrew S. Brown, Ph.D.',
    author_email='andrew@bridgeinformatics.com',
    description='RabidSeq Pipeline',
    url='www.violettx.com',
    install_requires=requirements,
    include_package_data=True
)
