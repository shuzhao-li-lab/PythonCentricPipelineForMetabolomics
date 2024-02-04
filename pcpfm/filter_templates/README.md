This directory has filters that can be passed to a variety of commands to exclude or include samples based on their metadata. 

For example 

```
{
    "Method": {
        "includes": [
            "HILIC_neg"
        ]
    }
}
```

Drops any sample that includes the substring "HILIC_neg" in the "Method" string. 


```
{
    "Method": {
        "excludes": [
            "HILIC_neg"
        ]
    }
}
```

While this will skip any sample that does not include the substring.