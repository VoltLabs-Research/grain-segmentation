cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(grain-segmentation VERSION 1.0.0 LANGUAGES CXX)

find_package(coretoolkit REQUIRED)
find_package(structure-identification REQUIRED)
find_package(polyhedral-template-matching REQUIRED)

volt_add_plugin(grain-segmentation
    DEPENDENCIES
        structure-identification::structure-identification
        polyhedral-template-matching::polyhedral-template-matching
)
