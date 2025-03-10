# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.0 - 2024-06-25
### Added
- add vidarr labels to outputs (changes to medata only).
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) 

## 1.0.7 - 2023-10-17
### Changed
- Update calculateContamination workflow to configure assembly-specific modules inside wdl.
- [GBS-4324](https://jira.oicr.on.ca/browse/GBS-4324) 

## 1.0.5 - 2022-10-13
### Added
- Added missing output prefix.

## 1.0.2 - 2022-05-31
### Changed
- Make tumor/normal sample selection explicit.
- Move hard-coded parameters to inputs.

## 1.0.1 - 2022-04-15
### Changed 
- Generating contamination metrics on T/N pairs
- [GP-3249](https://jira.oicr.on.ca/browse/GP-3249) 
- Accept (or create) bam/bai inputs to run `gatk calculateContamination` on.
- Designed for QC in other workflows

## 1.0.0 - 2022-04-27
### Changed 
- Some cleanups.