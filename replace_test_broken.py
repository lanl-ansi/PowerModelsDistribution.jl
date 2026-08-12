
from pathlib import Path

def find_and_replace(directory_path, search_text, replace_text, extension=None):
    # Convert string path to a Path object
    dir_path = Path(directory_path)
    
    # Define search pattern (all files or specific extension)
    search_pattern = f"*{extension}" if extension else "*"
    
    # rglob scans the directory and all subdirectories recursively
    for file_path in dir_path.rglob(search_pattern):
        if file_path.is_file() and str(file_path).endswith(".jl"):
            try:
                # Read the contents of the file
                content = file_path.read_text(encoding="utf-8")
                
                # Check if the target text exists to avoid unnecessary writes
                if search_text in content:
                    new_content = content.replace(search_text, replace_text)
                    
                    # Write the modified content back to the file
                    file_path.write_text(new_content, encoding="utf-8")
                    print(f"Updated: {file_path}")
                    
            except (UnicodeDecodeError, PermissionError):
                # Skip binary files or restricted files safely
                continue

# Example Usage:
find_and_replace("./test", "@test_broken", "@test")
# find_and_replace("/path/to/your/directory", "old_string", "new_string", extension=".txt")