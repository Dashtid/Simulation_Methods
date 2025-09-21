# Simulation Methods in Biomedical Engineering - Claude Code Guidelines

## Project Overview

This repository contains MATLAB code and scripts for various laboratory assignments in the course Simulation Methods in Biomedical Engineering. The project covers a range of numerical methods used in biomedical engineering, including solving partial differential equations (PDEs), boundary value problems, molecular dynamics simulations, and the finite element method (FEM). The implementations demonstrate practical applications of computational methods in biomedical engineering contexts.

## Development Environment

**Operating System**: Windows 11
**Shell**: Git Bash / PowerShell / Command Prompt
**Important**: Always use Windows-compatible commands:
- Use `dir` instead of `ls` for Command Prompt
- Use PowerShell commands when appropriate
- File paths use backslashes (`\`) in Windows
- Use `python -m http.server` for local development server
- Git Bash provides Unix-like commands but context should be Windows-aware

## Development Guidelines

### Code Quality
- Follow MATLAB coding conventions and best practices
- Use meaningful variable and function names
- Implement proper error handling for numerical computations
- Add comprehensive comments explaining algorithm steps and mathematical concepts
- Use consistent indentation and formatting
- Maintain clean, readable code
- Follow language-specific best practices

### Security
- No sensitive information in the codebase
- Use HTTPS for all external resources
- Regular dependency updates
- Follow security best practices for the specific technology stack

### MATLAB-Specific Guidelines for Numerical Methods
- Use vectorized operations where possible for performance
- Implement proper memory management for large computational problems
- Use appropriate MATLAB toolboxes (Partial Differential Equation, Optimization)
- Follow MATLAB naming conventions (camelCase for functions, lowercase for variables)
- Implement numerical stability checks and convergence criteria
- Use built-in MATLAB functions for linear algebra operations
- Include proper scaling and conditioning for numerical problems
- Implement appropriate visualization for simulation results

### Biomedical Engineering Specific Guidelines
- Document physical parameters and their units clearly
- Validate numerical results against analytical solutions where possible
- Include proper boundary condition implementations
- Use appropriate time stepping schemes for dynamic simulations
- Implement convergence studies for numerical accuracy
- Document assumptions and limitations of each simulation method

## Learning and Communication
- Always explain coding actions and decisions to help the user learn
- Describe why specific approaches or technologies are chosen
- Explain the purpose and functionality of code changes
- Provide context about best practices and coding patterns used
- Provide detailed explanations in the console when performing tasks, as many concepts may be new to the user