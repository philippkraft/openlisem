import os
import re

def generate_function_list(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    
    # Enhanced regular expression to capture functions belonging to the TWorld class
    function_regex = re.compile(
        r'[A-Za-z0-9*]*\s*[A-Za-z0-9_*]+\s*[A-Za-z0-9_~:*]+\(.*\)\s*\{\s*.*?\s*\}'
    )

    functions = []

    for line in lines:
        match = function_regex.match(line)
        if match:
            return_type, function_name, params, const = match.groups()
            function_signature = f"{return_type.strip()} {function_name}({params.strip()})"
            if const:
                function_signature += " const"
            functions.append(function_signature)
    
    if functions:
        functions_list = "/*\n * List of Functions:\n" + "\n".join(f" * {func};" for func in functions) + "\n */\n\n"
    else:
        functions_list = "/*\n * No functions found\n */\n\n"
    
    new_lines = [functions_list] + lines

    with open(file_path, 'w', encoding='utf-8') as file:
        file.writelines(new_lines)

def process_directory(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.cpp'):
                file_path = os.path.join(root, file)
                generate_function_list(file_path)


if __name__ == "__main__":
    #directory = input("Enter the directory path containing .cpp files: ")
    generate_function_list('c:/prgc/lisemgit/openlisem/test.cpp')
    #process_directory('c:/prgc/lisemgit/openlisem')
