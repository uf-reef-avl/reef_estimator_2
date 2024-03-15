cmake_minimum_required(VERSION 3.5)
project(reef_estimator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -I /usr/include/eigen3")
#############
# QTCREATOR #
#############
FILE(GLOB_RECURSE LibFiles "include/*")
add_custom_target(headers SOURCES ${LibFiles})

# Find dependencies
find_package(ament_cmake REQUIRED)
find_package(ament_cmake_python REQUIRED)
find_package(rclcpp REQUIRED)
find_package(rclpy REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(rosflight_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(std_msgs REQUIRED)
find_package(tf2_eigen REQUIRED)
find_package(tf_conversions REQUIRED)
find_package(reef_msgs REQUIRED)

# Include Cpp "include" directory
include_directories(include)

# Create Cpp executable
add_executable(
    reef_estimator
    src/reef_estimator_node.cpp
    src/sensor_manager.cpp
    src/estimator.cpp
    src/z_estimator.cpp
    src/xyz_estimator.cpp include/xyz_estimator.h
    src/xy_estimator.cpp
)
ament_target_dependencies(reef_estimator
  rclcpp
  std_msgs)

# Install Cpp executables
install(TARGETS reef_estimator
    DESTINATION lib/${PROJECT_NAME})
install(DIRECTORY include/
  DESTINATION include)
  
# Install Python modules
ament_python_install_package(${PROJECT_NAME})

# Install Python executables
install(PROGRAMS scripts/verify_estimates.py
        DESTINATION $lib/${PROJECT_NAME})

ament_export_dependencies(rclcpp rclpy sensor_msgs geometry_msgs rosflight rosflight_msgs sensor_msgs std_msgs tf2_eigen)
ament_export_include_directories()
ament_export_libraries()
ament_package()
