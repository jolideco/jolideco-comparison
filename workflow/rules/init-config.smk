rule init_config:
    input:
        "config/template.yaml"
    output:
        "config/{scenario}/{bkg_level}/{prefix}.yaml"
    log:
        "logs/init-config/{scenario}-{bkg_level}-{prefix}.log"
    run:
        import logging
        from pathlib import Path

        log = logging.getLogger(__name__)
        
        path = Path(input[0])
    
        with path.open("r") as f:
            template = f.read()

        max_value = MAX_VALUES[wildcards.scenario]
        title = f"{wildcards.scenario.title()} {wildcards.bkg_level.title()} {wildcards.prefix.title()}"

        template_filled = template.format(
                        title=title,
                        prefix=wildcards.prefix,
                        bkg_level=wildcards.bkg_level,
                        max_value=max_value,
                        scenario=wildcards.scenario
                    )

        path = Path(output[0])
       
        with path.open("w") as f:
            log.info(f"writing {path}")
            f.write(template_filled)