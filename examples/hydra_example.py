import hydra
from omegaconf import DictConfig

import logging

log = logging.getLogger(__name__)

@hydra.main(config_path="configs", config_name="default")
def some_test(cfg: DictConfig) -> None:

    log.info( f"{cfg.param1}, {cfg.param2}" )

if __name__ == "__main__":
    some_test()