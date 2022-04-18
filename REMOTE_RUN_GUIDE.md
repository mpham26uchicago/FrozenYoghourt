Guide to run module on remote machine
-------------------------------------

1. Create a jupyter notebook in the directory (same one as README).
2. Add this code at the beginning of the file and run.

```python
import sys
sys.path.insert(0, '')
from FrozenYoghourt import *

default_import()
```

3. The rest of the import code is in clipboard, simply pasted into the next cell